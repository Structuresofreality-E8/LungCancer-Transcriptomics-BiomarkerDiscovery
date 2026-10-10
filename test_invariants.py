"""Software invariants on synthetic inputs, never reported as cancer-study results."""
import unittest
import numpy as np
import pandas as pd
import luad_pipeline as p

class Invariants(unittest.TestCase):
    def test_patient_groups_never_cross_folds(self):
        groups=np.repeat(np.arange(30),2)
        y=np.tile([0,1],30);X=np.zeros((60,2))
        for train,test in p.splits(X,y,groups,5,42):
            self.assertFalse(set(groups[train]) & set(groups[test]))

    def test_sample_rank_is_scale_invariant_and_sample_local(self):
        x=pd.DataFrame([[0,1,5,5],[5,4,3,2]],columns=list('ABCD'))
        original=p.ranks(x)
        np.testing.assert_allclose(original,p.ranks(x*7))
        np.testing.assert_allclose(original.iloc[0],p.ranks(x.iloc[:1]).iloc[0])
        changed=x.copy();changed.iloc[1]=[100,0,0,0]
        np.testing.assert_allclose(original.iloc[0],p.ranks(changed).iloc[0])

    def test_variance_is_fit_on_training_only(self):
        train=np.array([[0,1],[0,2],[0,3]])
        fit=p.TrainingVariance().fit(train)
        test=np.array([[99,4],[101,5]])
        np.testing.assert_array_equal(fit.transform(test),[[4],[5]])

    def test_bh_and_missing_tests(self):
        np.testing.assert_allclose(p.bh([.01,.04,.03,np.nan]),[.03,.04,.04,np.nan],equal_nan=True)

if __name__=='__main__': unittest.main()
