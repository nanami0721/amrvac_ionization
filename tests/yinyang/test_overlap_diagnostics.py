"""Analytic checks for the common-physical-point diagnostic, no solver required."""
import unittest
import numpy as np
from overlap_diagnostics import basis, reconstruct


def fixture(n, patch, field):
    h = dict(xmin=np.array([1., .2*np.pi, -.8*np.pi]),
             xmax=np.array([2., .8*np.pi, .8*np.pi]), domain_nx=np.array([n, 2*n, 4*n]))
    axes = [a+(np.arange(k)+.5)*(b-a)/k for a,b,k in zip(h['xmin'],h['xmax'],h['domain_nx'])]
    r,t,p = np.meshgrid(*axes,indexing='ij')
    q = basis(t,p,patch)
    xyz = r[...,None]*q[..., :,0]
    w = np.zeros((*r.shape,8))
    w[...,5:8] = np.einsum('...ji,...j->...i',q,field(xyz))
    return h,w


class OverlapDiagnosticTests(unittest.TestCase):
    def test_uniform_global_vector_and_mask(self):
        field = lambda x: np.broadcast_to([.2,-.1,.3],x.shape)
        points=np.array([[1.,1.,.5],[-1.,1.,.5],[0.,0.,1.5],[4.,0.,0.]])
        for patch in (0,1):
            h,w=fixture(8,patch,field)
            b,valid=reconstruct(h,w,points,patch)
            self.assertTrue(valid.any())
            self.assertFalse(valid[-1])
            np.testing.assert_allclose(b[valid],field(points[valid]),atol=2e-15,rtol=0)
            self.assertTrue(np.isnan(b[~valid]).all())

    def test_smooth_reconstruction_converges(self):
        field=lambda x: np.stack((x[...,1],-x[...,0],np.full(x.shape[:-1],.3)),-1)
        points=np.array([[.8,.8,.8],[.9,.7,.9],[.8,.7,.9]])
        for patch in (0,1):
            errors=[]
            for n in (8,16,32):
                h,w=fixture(n,patch,field)
                b,valid=reconstruct(h,w,points,patch)
                self.assertTrue(valid.all())
                errors.append(float(np.linalg.norm(b-field(points),axis=1).mean()))
            self.assertLess(errors[1],errors[0]*.6)
            self.assertLess(errors[2],errors[1]*.6)


if __name__=='__main__':unittest.main()
