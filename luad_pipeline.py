"""LUAD tissue biomarker research pipeline. Run in Colab CPU, Python >=3.11.

Install first (fresh runtime):
%pip install -q pydeseq2==0.5.2 anndata==0.11.4 zarr==2.18.7 numcodecs==0.15.1 scikit-learn==1.6.1 gseapy==1.1.9 streamlit==1.45.1 pandas==2.2.3 numpy==2.2.6 scipy==1.15.3 matplotlib==3.10.3 seaborn==0.13.2 joblib==1.4.2
Then paste this entire file into the next cell, or use %run luad_pipeline.py.
Default run downloads the full eligible TCGA cohort; cached downloads resume.
"""
from pathlib import Path
import os, sys, json, re, io, gzip, hashlib, time, warnings, shutil
import urllib.request, urllib.parse, urllib.error
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import importlib.metadata
import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon, spearmanr, hypergeom
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, GridSearchCV
from sklearn.metrics import (roc_auc_score, average_precision_score, balanced_accuracy_score,
    brier_score_loss, confusion_matrix, RocCurveDisplay, PrecisionRecallDisplay)
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import joblib

# EDIT ONLY THIS CONFIGURATION. A persistent Drive path is useful for restarts.
ROOT = Path(os.environ.get('LUAD_ROOT', '/content/LungCancer-Transcriptomics-BiomarkerDiscovery'))
SEED = 42
DOWNLOAD_WORKERS = 3
N_JOBS = 2
OUTER_FOLDS = 5
INNER_FOLDS = 3
BOOTSTRAPS = 1000
RUN_DE = True
RUN_PATHWAYS = True
RUN_OPEN_TARGETS = True
REACTOME_LIBRARY = 'Reactome_2022'  # fixed database version; saved locally
MIN_COUNTS = 10
MIN_SAMPLES = 10
FEATURE_COUNTS = [10, 25, 50]  # upper bounds, not guaranteed signature sizes
GRID = {'select__k': FEATURE_COUNTS, 'clf__C': [0.01, 0.1, 1.0],
        'clf__l1_ratio': [0.5, 0.9]}
GEO_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40419/suppl/GSE40419_LC-87_RPKM_expression.txt.gz'
GEO_SOFT = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40419/soft/GSE40419_family.soft.gz'
API = 'https://api.gdc.cancer.gov'

def fetch(url, dest, md5=None):
    """Atomic, resumable-by-file download, retries, and GDC checksum checking."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and (not md5 or hashlib.md5(dest.read_bytes()).hexdigest() == md5):
        return dest
    for attempt in range(5):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'LUAD-research/1.0'})
            with urllib.request.urlopen(request, timeout=180) as src, open(str(dest)+'.part','wb') as dst:
                shutil.copyfileobj(src, dst)
            tmp = Path(str(dest)+'.part')
            if md5 and hashlib.md5(tmp.read_bytes()).hexdigest() != md5:
                raise ValueError('Checksum mismatch: '+str(dest))
            tmp.replace(dest)
            return dest
        except (OSError, ValueError, urllib.error.URLError):
            if attempt == 4:
                raise
            time.sleep(2**attempt)

def read_json(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except (OSError, urllib.error.URLError):
            if attempt == 4:
                raise
            time.sleep(2**attempt)

def bh(p):
    """BH over finite p-values; preserve unavailable tests as NaN."""
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if len(ok):
        order = ok[np.argsort(p[ok])]
        adj = p[order]*len(order)/np.arange(1,len(order)+1)
        q[order] = np.minimum(1, np.minimum.accumulate(adj[::-1])[::-1])
    return q

def save_table(df, name):
    df.to_csv(ROOT/'results'/name, index=True)

def figure(name):
    plt.tight_layout()
    plt.savefig(ROOT/'figures'/f'{name}.png', dpi=180, bbox_inches='tight')
    plt.savefig(ROOT/'figures'/f'{name}.pdf', bbox_inches='tight')
    plt.close()

def gdc_manifest():
    path = ROOT/'data/gdc_manifest.json'
    if path.exists():
        payload = json.loads(path.read_text())
    else:
        filt = {'op':'and','content':[
            {'op':'=','content':{'field':'cases.project.project_id','value':'TCGA-LUAD'}},
            {'op':'=','content':{'field':'data_type','value':'Gene Expression Quantification'}},
            {'op':'=','content':{'field':'analysis.workflow_type','value':'STAR - Counts'}},
            {'op':'=','content':{'field':'access','value':'open'}}]}
        params = {'filters':json.dumps(filt), 'format':'JSON', 'size':2000,
            'expand':'cases.samples,cases.diagnoses,cases.demographic',
            'fields':','.join(['file_id','file_name','md5sum','file_size','cases.submitter_id',
            'cases.case_id','cases.samples.sample_type','cases.samples.submitter_id',
            'cases.samples.sample_id','cases.diagnoses.ajcc_pathologic_stage',
            'cases.demographic.sex_at_birth','cases.demographic.race',
            'cases.demographic.age_at_index'])}
        payload = read_json(API+'/files?'+urllib.parse.urlencode(params))
        assert payload['data']['pagination']['total'] <= 2000, 'Implement pagination before proceeding.'
        path.write_text(json.dumps(payload,indent=2))
    rows, rejected = [], []
    for hit in payload['data']['hits']:
        cases = hit.get('cases',[])
        if len(cases)!=1 or len(cases[0].get('samples',[]))!=1:
            rejected.append({'file_id':hit['file_id'],'reason':'ambiguous case/sample association'})
            continue
        case, sample = cases[0], cases[0]['samples'][0]
        kind = sample.get('sample_type')
        if kind not in ['Primary Tumor','Solid Tissue Normal']:
            continue
        primary = [d for d in case.get('diagnoses',[]) if d.get('diagnosis_is_primary_disease') is not False]
        stages = sorted({d['ajcc_pathologic_stage'] for d in primary if d.get('ajcc_pathologic_stage')})
        demo = case.get('demographic',{})
        patient = case['submitter_id']
        rows.append({'file_id':hit['file_id'],'file_name':hit['file_name'],
            'md5sum':hit['md5sum'],'patient':patient,'sample':sample['submitter_id'],
            'condition':'Tumor' if kind=='Primary Tumor' else 'Normal',
            'stage':stages[0] if len(stages)==1 else 'unknown',
            'sex':demo.get('sex_at_birth',demo.get('gender','unknown')),
            'age':demo.get('age_at_index',np.nan),'race':demo.get('race','unknown'),
            'tss':patient.split('-')[1]})
    meta = pd.DataFrame(rows).sort_values(['patient','condition','sample','file_id'])
    dup = meta.duplicated(['patient','condition'],keep='first')
    # Lexical selection is label-independent, reproducible, and does not cherry-pick library size.
    save_table(meta.loc[dup], 'excluded_duplicate_aliquots.csv')
    save_table(pd.DataFrame(rejected), 'ambiguous_files.csv')
    meta = meta.loc[~dup].set_index('file_id')
    assert len(meta)>100 and (meta.condition=='Normal').sum()>=20, 'Unexpected TCGA cohort.'
    save_table(meta,'tcga_metadata.csv')
    save_table(meta.isna().sum().rename('null_count').to_frame(),'tcga_missing_metadata.csv')
    return meta

def read_star(path):
    frame = pd.read_csv(path,sep='\t',comment='#')
    required = {'gene_id','gene_name','gene_type','unstranded','fpkm_unstranded'}
    assert required <= set(frame.columns), f'Unexpected GDC columns: {frame.columns.tolist()}'
    frame = frame[frame.gene_id.str.startswith('ENSG')].copy()
    # Keep PAR_Y suffixes: naive split('.') merges X and pseudoautosomal Y loci.
    frame['gene_id'] = frame.gene_id.str.replace(r'\.\d+', '', regex=True)
    assert not frame.gene_id.duplicated().any(), 'Ambiguous Ensembl IDs.'
    return frame.set_index('gene_id')

def load_tcga(meta):
    def get(fid):
        return fetch(API+'/data/'+fid, ROOT/'data/gdc'/f'{fid}.tsv',meta.loc[fid,'md5sum'])
    ids = meta.index.tolist()
    paths = []
    with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as pool:
        for i,path in enumerate(pool.map(get,ids),1):
            paths.append(path)
            if i%25==0 or i==len(ids): print(f'TCGA downloads checked: {i}/{len(ids)}',flush=True)
    counts, fpkm, ann, metrics = {}, {}, None, []
    for fid, path in zip(ids,paths):
        frame = read_star(path)
        if ann is None:
            ann = frame[['gene_name','gene_type']].copy()
        else:
            assert frame.index.equals(ann.index), 'Annotation version mismatch; harmonize explicitly.'
            assert frame[['gene_name','gene_type']].equals(ann), 'Inconsistent GDC gene annotation.'
        c = frame.unstranded.to_numpy(float)
        assert np.isfinite(c).all() and (c>=0).all() and np.equal(c,np.floor(c)).all()
        counts[fid] = frame.unstranded.astype(np.int64)
        fpkm[fid] = frame.fpkm_unstranded.astype(float)
        raw = pd.read_csv(path,sep='\t',comment='#').set_index('gene_id')
        q = {'file_id':fid,'assigned_gene_reads':int(c.sum()), 'detected_genes':int((c>0).sum())}
        for name in ['N_unmapped','N_multimapping','N_noFeature','N_ambiguous']:
            q[name] = float(raw.loc[name,'unstranded']) if name in raw.index else np.nan
        metrics.append(q)
    counts = pd.DataFrame(counts).T.loc[ids]
    fpkm = pd.DataFrame(fpkm).T.loc[ids]
    assert np.isfinite(fpkm.values).all() and (fpkm.values>=0).all()
    qc = pd.DataFrame(metrics).set_index('file_id').join(meta)
    med = np.median(np.log1p(qc.assigned_gene_reads))
    mad = np.median(abs(np.log1p(qc.assigned_gene_reads)-med))
    qc['library_size_flag'] = abs(np.log1p(qc.assigned_gene_reads)-med)>4*max(1.4826*mad,0.01)
    # Flags are reported; no post-hoc deletion to improve classification.
    save_table(qc,'sequencing_summary_qc.csv')
    save_table(ann,'gdc_gene_annotation.csv')
    return counts, fpkm, ann, qc

def load_geo():
    path = fetch(GEO_URL,ROOT/'data/GSE40419_RPKM.txt.gz')
    soft = fetch(GEO_SOFT,ROOT/'data/GSE40419_family.soft.gz')
    # Use GEO metadata, not column ordering or assumed tumor-first layout.
    metadata, rec = [], None
    with gzip.open(soft,'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE = '):
                if rec is not None: metadata.append(rec)
                rec = {'gsm':line.split(' = ',1)[1],'characteristics':[]}
            elif rec is not None and line.startswith('!Sample_title = '):
                rec['title'] = line.split(' = ',1)[1]
            elif rec is not None and line.startswith('!Sample_characteristics_ch1 = '):
                rec['characteristics'].append(line.split(' = ',1)[1])
        if rec is not None: metadata.append(rec)
    gm = pd.DataFrame(metadata)
    def label(title):
        # Actual GEO schema: LC_C1, LC_S1, LC_C1_nor, LC_S1_nor.
        m = re.fullmatch(r'LC_([CS]\d+)(_nor)?',str(title),flags=re.I)
        assert m, f'Unrecognized GEO sample: {title}. Inspect SOFT before changing parser.'
        return 'GEO_'+m[1].upper(), 'Normal' if m[2] else 'Tumor'
    gm[['patient','condition']] = pd.DataFrame([label(t) for t in gm.title],index=gm.index)
    gm = gm.set_index('title')
    assert gm.index.is_unique, 'Duplicate GEO sample titles need explicit resolution.'
    df = pd.read_csv(path,sep='\t',low_memory=False)
    assert 'gene' in df.columns, f'Unexpected GEO columns: {df.columns[:10].tolist()}'
    sample_cols = [c for c in df.columns if c in gm.index]
    assert len(sample_cols)==164 and len(gm)==164, 'Expected all 164 GEO expression samples.'
    gm = gm.loc[sample_cols]
    assert gm.condition.value_counts().to_dict()=={'Tumor':87,'Normal':77}
    values = df[sample_cols].apply(pd.to_numeric,errors='raise')
    assert np.isfinite(values.values).all() and (values.values>=0).all()
    genes = df.gene.astype(str).str.strip()
    valid = df.gene.notna() & ~genes.isin(['','-','NA','nan'])
    # Transcript duplicates: median RPKM, fixed in advance, never sum transcript RPKM.
    save_table(pd.DataFrame({'symbol':genes[valid].value_counts()}),'geo_symbol_multiplicity.csv')
    values.index = genes
    expr = values.loc[valid.to_numpy()].groupby(level=0,sort=True).median().T
    assert not expr.columns.duplicated().any()
    save_table(gm,'geo_metadata.csv')
    return expr, gm

def paired_ids(meta):
    patients = meta.groupby('patient').condition.nunique()
    return meta.index[meta.patient.isin(patients[patients==2].index)]

def paired_de(counts,meta,ann):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    ids = paired_ids(meta)
    assert len(ids)>=20, 'Insufficient complete patient pairs.'
    c = counts.loc[ids]
    keep = (c>=MIN_COUNTS).sum(axis=0)>=MIN_SAMPLES
    # Protein-coding unique symbols; same annotation-only restriction used for modeling.
    unique = ann.gene_type.eq('protein_coding') & ann.gene_name.notna() & ~ann.gene_name.duplicated(keep=False)
    c = c.loc[:,keep & unique]
    md = meta.loc[ids,['patient','condition']].copy()
    md['condition'] = pd.Categorical(md.condition,categories=['Normal','Tumor'])
    dds = DeseqDataSet(counts=c,metadata=md,design='~patient+condition',
                       refit_cooks=True,n_cpus=N_JOBS,quiet=True,low_memory=True)
    dds.deseq2()
    stats = DeseqStats(dds,contrast=['condition','Tumor','Normal'],alpha=0.05,n_cpus=N_JOBS,quiet=True)
    stats.summary()
    result = stats.results_df.join(ann)
    save_table(result,'paired_differential_expression.csv')
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers['vst_counts'],index=c.index,columns=c.columns)
    pcs = PCA(n_components=2,random_state=SEED).fit_transform(vst)
    plt.figure(figsize=(6,5));sns.scatterplot(x=pcs[:,0],y=pcs[:,1],hue=md.condition)
    plt.xlabel('PC1');plt.ylabel('PC2');plt.title('Paired TCGA: blind VST PCA')
    figure('paired_vst_pca')
    plt.figure(figsize=(7,5))
    sns.scatterplot(x=result.log2FoldChange,y=-np.log10(result.padj.clip(lower=1e-300)),
                    hue=result.padj.lt(0.05),s=8,linewidth=0)
    plt.xlabel('Tumor / normal log2 fold change');plt.ylabel('-log10 BH-adjusted p')
    plt.title('Paired differential expression');figure('volcano')
    return result

def ranks(expr):
    # Per-sample percentile ranks over a fixed gene universe: no cohort-wide fitting.
    # FPKM and RPKM both account for gene length; ranks reduce scale shifts, not all protocol bias.
    return pd.DataFrame(rankdata(expr.to_numpy(float),axis=1,method='average')/expr.shape[1],
                        index=expr.index,columns=expr.columns)

class TrainingVariance(BaseEstimator, TransformerMixin):
    def fit(self,X,y=None):
        self.keep_ = np.var(X,axis=0)>1e-12
        if not self.keep_.any(): raise ValueError('No variable genes in training split.')
        return self
    def transform(self,X): return np.asarray(X)[:,self.keep_]

class StrictLogisticRegression(LogisticRegression):
    def fit(self,X,y,sample_weight=None):
        from sklearn.exceptions import ConvergenceWarning
        with warnings.catch_warnings():
            warnings.simplefilter('error',ConvergenceWarning)
            return super().fit(X,y,sample_weight=sample_weight)

def model():
    return Pipeline([('variance',TrainingVariance()),('select',SelectKBest(f_classif)),
        ('scale',StandardScaler()),('clf',StrictLogisticRegression(penalty='elasticnet',solver='saga',
            class_weight='balanced',max_iter=15000,tol=1e-4,random_state=SEED))])

def splits(X,y,g,n,seed):
    folds = list(StratifiedGroupKFold(n_splits=n,shuffle=True,random_state=seed).split(X,y,g))
    for tr,te in folds:
        assert set(g[tr]).isdisjoint(set(g[te])), 'Patient leakage.'
        assert len(np.unique(y[tr]))==len(np.unique(y[te]))==2, 'Single-class split; reduce folds.'
    return folds

def search(X,y,g,seed=SEED):
    cv = splits(X,y,g,INNER_FOLDS,seed)
    grid = GridSearchCV(model(),GRID,scoring='roc_auc',cv=cv,n_jobs=N_JOBS,error_score='raise',refit=True)
    grid.fit(X,y)
    return grid

def coefficients(pipe,genes):
    available = np.asarray(genes)[pipe.named_steps['variance'].keep_]
    selected = available[pipe.named_steps['select'].get_support()]
    coeff = pd.Series(pipe.named_steps['clf'].coef_[0],index=selected)
    return coeff[abs(coeff)>1e-8]

def metrics(y,p):
    pred = p>=0.5
    tn,fp,fn,tp = confusion_matrix(y,pred,labels=[0,1]).ravel()
    return {'roc_auc':float(roc_auc_score(y,p)), 'average_precision':float(average_precision_score(y,p)),
        'balanced_accuracy':float(balanced_accuracy_score(y,pred)),
        'sensitivity':float(tp/(tp+fn)), 'specificity':float(tn/(tn+fp)),
        'brier_uncalibrated':float(brier_score_loss(y,p)), 'tumor_prevalence':float(np.mean(y)),
        'n':int(len(y)), 'tn':int(tn),'fp':int(fp),'fn':int(fn),'tp':int(tp)}

def patient_ci(y,p,groups):
    rng = np.random.default_rng(SEED)
    unique = np.unique(groups)
    blocks = {g:np.flatnonzero(groups==g) for g in unique}
    samples = []
    for _ in range(BOOTSTRAPS):
        idx = np.concatenate([blocks[g] for g in rng.choice(unique,len(unique),replace=True)])
        if len(np.unique(y[idx]))==2: samples.append(metrics(y[idx],p[idx]))
    return {key:[float(x) for x in np.quantile([s[key] for s in samples],[.025,.975])]
            for key in ['roc_auc','average_precision','balanced_accuracy','sensitivity','specificity']}

def nested_fit(expr,meta,qc):
    X = expr.values
    y = meta.condition.eq('Tumor').to_numpy(int)
    groups = meta.patient.to_numpy()
    outer = splits(X,y,groups,OUTER_FOLDS,SEED)
    p = np.full(len(y),np.nan);baseline = p.copy();fold_id = np.zeros(len(y),int)
    coefs,params = [],[]
    qx = np.log1p(qc.loc[meta.index,['assigned_gene_reads','detected_genes']].values)
    for i,(tr,te) in enumerate(outer):
        print(f'Nested CV fold {i+1}/{OUTER_FOLDS}',flush=True)
        fit = search(X[tr],y[tr],groups[tr],SEED+i)
        p[te] = fit.predict_proba(X[te])[:,1];fold_id[te]=i+1
        coefs.append(coefficients(fit.best_estimator_,expr.columns))
        params.append({'fold':i+1,**fit.best_params_,'test_auc':roc_auc_score(y[te],p[te])})
        b = Pipeline([('s',StandardScaler()),('m',LogisticRegression(class_weight='balanced'))])
        b.fit(qx[tr],y[tr]);baseline[te] = b.predict_proba(qx[te])[:,1]
    assert np.isfinite(p).all()
    save_table(pd.DataFrame(params),'outer_fold_parameters.csv')
    save_table(meta.assign(tumor_score=p,fold=fold_id,qc_only_score=baseline),'tcga_oof_predictions.csv')
    stability = pd.concat(coefs,axis=1).reindex(expr.columns).fillna(0)
    stability.columns = [f'fold_{i+1}' for i in range(OUTER_FOLDS)]
    stability['selection_frequency'] = (abs(stability)>1e-8).mean(axis=1)
    save_table(stability,'outer_fold_gene_stability.csv')
    final = search(X,y,groups)
    sig = coefficients(final.best_estimator_,expr.columns).rename('standardized_coefficient').to_frame()
    sig = sig.join(stability[['selection_frequency']]).sort_values('standardized_coefficient')
    assert len(sig)>0, 'Regularization produced an intercept-only model; report failure, do not inspect GEO to fix it.'
    save_table(sig,'locked_signature.csv')
    # Lock model and feature universe BEFORE accessing GEO labels for evaluation.
    bundle = {'pipeline':final.best_estimator_,'genes':expr.columns.tolist(),
        'normalization':'within-sample percentile rank over all saved genes, average ties',
        'best_params':final.best_params_,'threshold':0.5,'signature':sig}
    joblib.dump(bundle,ROOT/'results/locked_model.joblib')
    pd.Series(expr.columns,name='gene').to_csv(ROOT/'results/fixed_gene_universe.csv',index=False)
    return bundle,p,baseline,sig

def replication(geo,gm,sig,de):
    ids = paired_ids(gm)
    md = gm.loc[ids]
    rows = []
    for gene in sig.index:
        tmp = pd.DataFrame({'value':np.log2(geo.loc[ids,gene]+1),
            'patient':md.patient,'condition':md.condition}).pivot(index='patient',columns='condition',values='value')
        delta = tmp.Tumor-tmp.Normal
        p = wilcoxon(delta,zero_method='wilcox').pvalue if (delta!=0).any() else 1.0
        row = {'gene':gene,'n_pairs':len(delta),'median_paired_log2_RPKM_plus1_difference':delta.median(),
               'paired_wilcoxon_p':p}
        rows.append(row)
    rep = pd.DataFrame(rows).set_index('gene')
    rep['paired_wilcoxon_BH'] = bh(rep.paired_wilcoxon_p)
    if de is not None:
        lookup = de.set_index('gene_name')[['log2FoldChange','padj']]
        rep = rep.join(lookup.rename(columns={'padj':'tcga_paired_DE_BH'}))
        rep['same_direction'] = np.sign(rep.log2FoldChange)==np.sign(rep.median_paired_log2_RPKM_plus1_difference)
        rep['replicated_signature_association'] = rep.same_direction & rep.tcga_paired_DE_BH.lt(.05) & rep.paired_wilcoxon_BH.lt(.05)
    save_table(rep,'external_paired_gene_replication.csv')
    return rep

def subgroup_report(meta,y,p,name):
    rows = []
    for field in ['sex','race','tss','stage']:
        if field not in meta: continue
        for val,idx in meta.groupby(field).indices.items():
            idx = np.asarray(idx)
            if len(idx)>=20 and min(np.bincount(y[idx],minlength=2))>=5:
                rows.append({'field':field,'group':val,**metrics(y[idx],p[idx])})
    save_table(pd.DataFrame(rows),name+'_subgroup_metrics.csv')

def composition_diagnostics(expr,meta,p,name):
    panels = {'immune':['PTPRC','CD3D','CD3E','MS4A1','LST1'],
              'stromal':['COL1A1','COL1A2','DCN','LUM'],
              'epithelial':['EPCAM','KRT8','KRT18','KRT19']}
    rows = []
    for panel,markers in panels.items():
        present = [g for g in markers if g in expr]
        if len(present)<2: continue
        proxy = expr[present].mean(axis=1)
        for label in ['Normal','Tumor']:
            take = meta.condition.eq(label).values
            if take.sum()<10: continue
            rho,pv = spearmanr(proxy.values[take],p[take])
            rows.append({'panel':panel,'condition':label,'markers':','.join(present),
                         'n':int(take.sum()),'rho_score_vs_proxy':rho,'p':pv})
    tab = pd.DataFrame(rows)
    if len(tab): tab['BH'] = bh(tab.p)
    save_table(tab,name+'_composition_proxy_diagnostics.csv')

def plots(tcga,tm,geo,gm,p,pe,sig,rep):
    for name,ex,md in [('tcga',tcga,tm),('geo',geo,gm)]:
        pcs = PCA(n_components=2,random_state=SEED).fit_transform(ex.values)
        plt.figure(figsize=(6,5));sns.scatterplot(x=pcs[:,0],y=pcs[:,1],hue=md.condition)
        plt.title(name.upper()+' sample rank PCA');plt.xlabel('PC1');plt.ylabel('PC2');figure(name+'_rank_pca')
        top = sig.standardized_coefficient.abs().nlargest(min(20,len(sig))).index
        order = md.sort_values(['condition','patient']).index
        # Training-derived scaling only, including the external heatmap.
        z = (ex.loc[order,top]-tcga[top].mean())/tcga[top].std().replace(0,1)
        colors = md.loc[order,'condition'].map({'Tumor':'#d95f02','Normal':'#1b9e77'})
        cg = sns.clustermap(z.T,row_cluster=True if len(top)>1 else False,col_cluster=False,
            col_colors=colors,vmin=-3,vmax=3,cmap='vlag',xticklabels=False,figsize=(11,6))
        cg.fig.suptitle(name.upper()+' locked signature: green normal, orange tumor',y=1.02)
        cg.savefig(ROOT/'figures'/f'{name}_signature_heatmap.png',dpi=180,bbox_inches='tight');plt.close(cg.fig)
    fig,axes = plt.subplots(1,2,figsize=(11,4))
    for label,y,pr in [('TCGA nested OOF',tm.condition.eq('Tumor'),p),('GEO independent',gm.condition.eq('Tumor'),pe)]:
        RocCurveDisplay.from_predictions(y,pr,name=label,ax=axes[0])
        PrecisionRecallDisplay.from_predictions(y,pr,name=label,ax=axes[1])
    figure('validation_curves')
    plt.figure(figsize=(7,max(3,len(sig)*.22)))
    plt.barh(sig.index,sig.standardized_coefficient);plt.axvline(0,color='black',lw=.7)
    plt.xlabel('Locked TCGA standardized coefficient; no GEO refitting');figure('locked_coefficients')
    plt.figure(figsize=(6,5))
    if 'log2FoldChange' in rep:
        sns.scatterplot(data=rep,x='log2FoldChange',y='median_paired_log2_RPKM_plus1_difference')
        plt.axhline(0,color='grey');plt.axvline(0,color='grey')
        plt.title('Signature direction replication; effect scales differ');figure('paired_gene_replication')
    else: plt.close()

def pathways(de):
    import gseapy as gp
    gmt = ROOT/'data'/f'{REACTOME_LIBRARY}_sets.json'
    if gmt.exists(): sets = json.loads(gmt.read_text())
    else:
        sets = gp.get_library(name=REACTOME_LIBRARY,organism='Human')
        gmt.write_text(json.dumps(sets,indent=2))
    rank = de.dropna(subset=['stat']).set_index('gene_name')['stat'].sort_values(ascending=False)
    assert rank.index.is_unique
    result = gp.prerank(rnk=rank,gene_sets=sets,min_size=15,max_size=500,permutation_num=1000,
                       threads=N_JOBS,seed=SEED,outdir=str(ROOT/'results/gsea'),no_plot=True)
    tab = result.res2d.copy();save_table(tab,'reactome_preranked_GSEA.csv')
    # Local ORA uses all tested finite-p genes as the background, not the entire genome.
    background = set(de.loc[de.pvalue.notna(),'gene_name'])
    significant = set(de.loc[de.padj.lt(.05)&de.log2FoldChange.abs().ge(1),'gene_name']) & background
    rows = []
    for term,genes in sets.items():
        members = set(genes)&background
        if not 15<=len(members)<=500: continue
        overlap = significant & members
        pv = hypergeom.sf(len(overlap)-1,len(background),len(members),len(significant))
        rows.append({'pathway':term,'p':pv,'overlap':len(overlap),'tested_members':len(members),
                     'genes':';'.join(sorted(overlap))})
    ora = pd.DataFrame(rows);ora['BH'] = bh(ora.p);ora = ora.sort_values('BH')
    save_table(ora,'reactome_ORA.csv')
    top = tab.assign(FDR=pd.to_numeric(tab['FDR q-val']),NES=pd.to_numeric(tab['NES'])).nsmallest(15,'FDR')
    plt.figure(figsize=(9,6));sns.scatterplot(data=top,x='NES',y='Term',size=-np.log10(top.FDR.clip(lower=1e-6)),hue='NES',palette='coolwarm')
    plt.title('Reactome preranked GSEA (gene-set permutation)');figure('reactome_dotplot')

def open_targets(sig,ann):
    # Post-validation annotation only. A single batched query for this small signature.
    unique = ann.loc[~ann.gene_name.duplicated(keep=False)].reset_index().set_index('gene_name')
    ids = [unique.loc[g,'gene_id'] for g in sig.index]
    fields = []
    for i,gid in enumerate(ids):
        fields.append('t'+str(i)+': target(ensemblId: '+json.dumps(gid)+') { id approvedSymbol '
            'associatedDiseases(page: {index: 0, size: 200}) { count rows { score disease { id name } } } }')
    payload = {'query':'{ '+' '.join(fields)+' }'}
    cache = ROOT/'data/open_targets_signature.json'
    cached = json.loads(cache.read_text()) if cache.exists() else {}
    if cached.get('query')==payload['query']:
        response = cached['response']
    else:
        req = urllib.request.Request('https://api.platform.opentargets.org/api/v4/graphql',
              data=json.dumps(payload).encode(),headers={'Content-Type':'application/json'})
        with urllib.request.urlopen(req,timeout=90) as r: response = json.load(r)
        if response.get('errors'): raise ValueError(response['errors'])
        cache.write_text(json.dumps({'query':payload['query'],'response':response,
            'queried_utc':datetime.now(timezone.utc).isoformat()},indent=2))
    rows = []
    for i,g in enumerate(sig.index):
        target = response['data'].get('t'+str(i))
        if target is None:
            rows.append({'gene':g,'status':'target not available'});continue
        assoc = target['associatedDiseases']
        lung = [r for r in assoc['rows'] if 'lung' in r['disease']['name'].lower()]
        if not lung:
            rows.append({'gene':g,'target_id':ids[i],
                'status':'no lung term in first 200 associations; NOT evidence of absence',
                'association_count':assoc['count']})
        for r in lung:
            rows.append({'gene':g,'target_id':ids[i],'disease_id':r['disease']['id'],
                'disease_name':r['disease']['name'],'association_score':r['score'],
                'association_count':assoc['count'],'status':'contextual annotation only'})
    save_table(pd.DataFrame(rows),'open_targets_lung_associations.csv')

LIMITATIONS = '''# Interpretation and next evidence required

This experiment classifies resected LUAD tissue versus adjacent nonmalignant
lung tissue. Adjacent normal is not a healthy screening population. Tissue
discrimination does not establish blood detectability or preclinical detection.

Mitigated: patient leakage (all folds grouped), global supervised feature
selection (fold-local), cross-study units (per-sample gene ranks), duplicate
aliquots (fixed lexical choice), DE pairing (complete-pair design), multiplicity
(BH), and external test tuning (locked model). Ranks do not remove all platform,
annotation, gene-length, cellular-composition, or population differences.

Patient bootstrap intervals describe test-patient sampling uncertainty at a
fixed trained model; nested OOF intervals do not capture full training-process
uncertainty. Average precision depends on prevalence. Class weighting means
scores and Brier values are not calibrated clinical probabilities. Threshold
0.5 was fixed before external testing. Bootstrap all model training stages or
repeat nested CV for a stronger, more expensive uncertainty analysis.

Signature stability is a descriptive frequency across outer training folds;
it is not a probability that a gene is causal. Replication p-values are BH
adjusted over the locked selected signature only. External log2(RPKM+1)
paired differences are not directly comparable to DESeq2 log2 fold changes.
GSEA uses gene-set permutations; gene correlation can affect significance.

Bulk immune/stromal/epithelial marker averages are confounding diagnostics,
not estimated cell fractions or purity. Confirm using pathology-derived purity,
validated deconvolution, and independent single-cell or spatial data. Do not
interpret signature genes as tumor-cell-intrinsic or causal without that work.

Age, smoking, library preparation, and medication may be incompletely captured.
GDC tissue source site is a center proxy, not a sequencing batch. Sparse
subgroup estimates are omitted, not proof of fairness. Basic STAR summary
metrics are not FASTQ QC; raw-read processing is required for read-level checks.

Next study: lock the signature, assay and threshold; test blinded in additional
independent tissue cohorts including benign nodules, inflammatory disease,
other cancers and metastatic lesions. For a blood-screening claim, assay the
markers in blood, assess detection limits and tissue-to-blood transport, then
validate in a prospective intended-use screening cohort with healthy controls,
benign disease, appropriate prevalence and stage distribution. Measure
sensitivity at prespecified high specificity, predictive values, calibration,
and clinical benefit. Neither TCGA nor GSE40419 can supply that evidence.

Stage-I tissue sensitivity is descriptive and does not demonstrate early
detection. Pathway/disease annotation is contextual evidence, not clinical
validation or novelty. No causal, diagnostic, prognostic, or treatment claim.
'''

def write_dashboard():
    app = '''from pathlib import Path
import streamlit as st
import pandas as pd
import json
ROOT = Path(__file__).resolve().parent
st.set_page_config(page_title="LUAD research explorer",layout="wide")
st.title("LUAD tissue biomarker research explorer")
st.caption("Retrospective research results; scores are not clinical probabilities.")
st.json(json.loads((ROOT/"results/metrics.json").read_text()))
sig = pd.read_csv(ROOT/"results/locked_signature.csv",index_col=0)
st.subheader("Locked signature")
st.dataframe(sig)
rep = pd.read_csv(ROOT/"results/external_paired_gene_replication.csv",index_col=0)
gene = st.selectbox("Inspect gene",sig.index.tolist())
st.dataframe(rep.loc[[gene]])
for name in ["validation_curves", "locked_coefficients", "geo_signature_heatmap", "reactome_dotplot"]:
    p = ROOT/"figures"/(name+".png")
    if p.exists(): st.image(str(p))
st.markdown((ROOT/"LIMITATIONS.md").read_text())
'''
    (ROOT/'app.py').write_text(app)

def main():
    for folder in ['data','results','figures']: (ROOT/folder).mkdir(parents=True,exist_ok=True)
    config = {k:globals()[k] for k in ['SEED','OUTER_FOLDS','INNER_FOLDS','BOOTSTRAPS','RUN_DE',
        'RUN_PATHWAYS','RUN_OPEN_TARGETS','REACTOME_LIBRARY','MIN_COUNTS','MIN_SAMPLES','GRID']}
    config['started_utc'] = datetime.now(timezone.utc).isoformat()
    (ROOT/'config.json').write_text(json.dumps(config,indent=2))
    versions = {p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','scikit-learn','matplotlib','joblib']}
    (ROOT/'runtime_versions.json').write_text(json.dumps(versions,indent=2))
    print('1/7 Downloading public expression files and metadata',flush=True)
    tm = gdc_manifest();counts,fpkm,ann,qc = load_tcga(tm)
    geo,gm = load_geo()
    # Gene availability only; no external values/labels drive filtering or tuning.
    unique = ann.gene_type.eq('protein_coding') & ann.gene_name.notna() & ~ann.gene_name.duplicated(keep=False)
    selected = ann.loc[unique & ann.gene_name.isin(geo.columns)]
    genes = sorted(selected.gene_name.tolist())
    assert len(genes)>5000, 'Insufficient unambiguous shared gene coverage; inspect mappings.'
    tc = fpkm.loc[:,selected.index].copy();tc.columns = selected.gene_name
    tc = tc[genes];ge = geo[genes]
    tc = ranks(tc);ge = ranks(ge)
    print(f'Shared fixed gene universe: {len(genes)}; TCGA {len(tm)}; GEO {len(gm)}',flush=True)
    print('2/7 Paired raw-count differential expression',flush=True)
    de = paired_de(counts,tm,ann) if RUN_DE else None
    print('3/7 Nested patient-grouped model fitting',flush=True)
    bundle,p,baseline,sig = nested_fit(tc,tm,qc)
    print('4/7 Evaluating locked model in GSE40419',flush=True)
    pe = bundle['pipeline'].predict_proba(ge.values)[:,1]
    yt = tm.condition.eq('Tumor').values.astype(int);ye = gm.condition.eq('Tumor').values.astype(int)
    report = {'tcga_nested_oof':metrics(yt,p),'tcga_oof_patient_bootstrap_CI':patient_ci(yt,p,tm.patient.values),
        'tcga_qc_only_baseline':metrics(yt,baseline),'geo_locked_external':metrics(ye,pe),
        'geo_patient_bootstrap_CI':patient_ci(ye,pe,gm.patient.values),
        'signature_size':len(sig),'shared_gene_count':len(genes),'score_threshold':0.5,
        'caution':'Tissue association only; class-weighted scores uncalibrated; CI conditional on fitted models.'}
    normal = tm.condition.eq('Normal')
    stage1 = tm.condition.eq('Tumor') & tm.stage.str.match(r'^stage i([abc])?$',case=False)
    take = (normal|stage1).values
    if stage1.sum()>=10:
        report['tcga_stage_I_tissue_vs_all_adjacent_normals_OOF'] = metrics(yt[take],p[take])
    (ROOT/'results/metrics.json').write_text(json.dumps(report,indent=2))
    save_table(gm.assign(tumor_score=pe,predicted_tumor=pe>=.5),'geo_locked_predictions.csv')
    rep = replication(geo[genes],gm,sig,de)
    subgroup_report(tm,yt,p,'tcga')
    composition_diagnostics(tc,tm,p,'tcga');composition_diagnostics(ge,gm,pe,'geo')
    print('5/7 Figures and pathway interpretation',flush=True)
    plots(tc,tm,ge,gm,p,pe,sig,rep)
    status = {'DE':'completed' if RUN_DE else 'disabled','pathways':'disabled','open_targets':'disabled'}
    if RUN_PATHWAYS and de is not None:
        try:
            pathways(de);status['pathways']='completed'
        except Exception as e:
            status['pathways']='FAILED: '+repr(e)
            warnings.warn('Pathway module failed; saved explicit status. No pathway conclusions are supported.')
    if RUN_OPEN_TARGETS:
        try:
            open_targets(sig,ann);status['open_targets']='completed'
        except Exception as e:
            status['open_targets']='FAILED: '+repr(e)
            warnings.warn('Open Targets annotation failed; saved explicit status.')
    (ROOT/'results/module_status.json').write_text(json.dumps(status,indent=2))
    print('6/7 Writing research explorer and reproducibility records',flush=True)
    (ROOT/'LIMITATIONS.md').write_text(LIMITATIONS)
    (ROOT/'.gitignore').write_text('data/gdc/\n*.zip\n__pycache__/\n')
    write_dashboard()
    # Save this source when run as a file; notebook stores the original source itself.
    source = globals().get('__file__')
    if source and Path(source).is_file() and Path(source).resolve()!= (ROOT/'luad_pipeline.py').resolve():
        shutil.copy2(source,ROOT/'luad_pipeline.py')
    all_versions = subprocess_versions()
    (ROOT/'requirements-lock.txt').write_text(all_versions)
    (ROOT/'README.md').write_text('''# LungCancer-Transcriptomics-BiomarkerDiscovery
Retrospective TCGA-LUAD to GSE40419 tissue association study.

Run in a fresh Colab CPU runtime using the supplied notebook. Default run uses
the full eligible cohort. Download progress is cached per file and checksum
verified. For persistence set LUAD_ROOT to a mounted Google Drive folder before
running. See config.json and runtime_versions.json for the actual run.

Model: shared unambiguous protein-coding gene symbols, length-normalized
expression, per-sample ranks, fold-local variance filtering and ANOVA selection,
scaling, class-weighted elastic-net logistic regression. Nested grouped 5x3 CV.
DE: raw unstranded counts on complete pairs with patient + condition design.
GEO duplicate transcript symbols: median RPKM. No external model tuning.
External values are used only after model fitting; public annotation-level
feature availability establishes the fixed measurement universe in advance.

Run explorer locally: `streamlit run app.py`. This is a results explorer,
not a clinical prediction application. No public tunnel is needed.

The archive includes cached public expression data for reproducibility. Do not
commit these bulk files to GitHub; publish source, manifests, configuration,
tables and figures. Review missing metadata, QC flags, convergence warnings,
module_status.json and LIMITATIONS.md before claiming results. Nonconvergence
requires rerunning with more iterations using TCGA only, before GEO evaluation.

Sources: https://portal.gdc.cancer.gov/ ;
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419 ;
https://pydeseq2.readthedocs.io/en/v0.5.2/ ;
Reactome via pinned Enrichr Reactome_2022 sets, cached in data/Reactome_2022_sets.json.
Open Targets lung-term associations are optional post-validation annotation
from the first 200 associations per target. They do not establish model validity
or novelty; review module_status.json to confirm API success.
''')
    print('7/7 Packaging outputs',flush=True)
    hashes = {str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
        for f in (ROOT/'data').rglob('*') if f.is_file()}
    (ROOT/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
    archive = shutil.make_archive(str(ROOT)+'_results','zip',ROOT)
    print(json.dumps(report,indent=2));print('Outputs:',ROOT);print('ZIP:',archive)
    return report

def subprocess_versions():
    import subprocess
    return subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True)

if __name__=='__main__':
    main()
