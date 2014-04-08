#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2011)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

import logging
import os
import os.path
from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
import sherpa.astro.sim as sim

from sherpa.astro.instrument import Response1D
from sherpa.astro.io import read_pha
from sherpa.astro.data import DataPHA
from sherpa.astro.xspec import XSwabs, XSpowerlaw
from sherpa.fit import Fit
from sherpa.stats import Cash, CStat
from sherpa.optmethods import NelderMead
from sherpa.estmethods import Covariance


logger = logging.getLogger('sherpa')


class test_sim(SherpaTestCase):

    def setUp(self):
        #self.startdir = os.getcwd()
        self.old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.CRITICAL)

        datadir = SherpaTestCase.datadir
        if datadir is None:
            return

        pha = os.path.join(datadir, "refake_0934_1_21_1e4.fak")
        rmf = os.path.join(datadir, "ccdid7_default.rmf")
        arf = os.path.join(datadir, "quiet_0934.arf")
        
        self.simarf = os.path.join(datadir, "aref_sample.fits")
        self.pcaarf = os.path.join(datadir, "aref_Cedge.fits")

        data = read_pha(pha)
        data.ignore(None,0.3)
        data.ignore(7.0,None)

        rsp = Response1D(data)
        self.abs1 = XSwabs('abs1')
        self.p1 = XSpowerlaw('p1')
        model = rsp(self.abs1*self.p1)
        
        self.fit = Fit(data, model, CStat(), NelderMead(), Covariance())


    def tearDown(self):
        #os.chdir(self.startdir)
        logger.setLevel(self.old_level)


    @needs_data
    def test_pragbayes_simarf(self):
        datadir = SherpaTestCase.datadir
        if datadir is None:
            return

        mcmc = sim.MCMC()

        self.abs1.nh = 0.092886
        self.p1.phoindex = 0.994544
        self.p1.norm = 9.26369

        mcmc.set_sampler("PragBayes")
        mcmc.set_sampler_opt("simarf", self.simarf)
        mcmc.set_sampler_opt("p_M", 0.5)
        mcmc.set_sampler_opt("nsubiter", 7)

        covar_results = self.fit.est_errors()
        cov = covar_results.extra_output

        niter = 10
        stats, accept, params = mcmc.get_draws(self.fit, cov, niter=niter)
        # try:
        #     assert (covar_results.parmaxes < params.std(1)).all()
        # except AssertionError:
        #     print 'covar: ', str(covar_results.parmaxes)
        #     print 'param: ', str(params.std(1))
        #     raise


    @needs_data
    def test_pragbayes_pcaarf(self):
        datadir = SherpaTestCase.datadir
        if datadir is None:
            return

        mcmc = sim.MCMC()

        self.abs1.nh = 0.092886
        self.p1.phoindex = 0.994544
        self.p1.norm = 9.26369

        mcmc.set_sampler("pragBayes")
        mcmc.set_sampler_opt("simarf", self.pcaarf)
        mcmc.set_sampler_opt("p_M", 0.5)
        mcmc.set_sampler_opt("nsubiter", 5)

        covar_results = self.fit.est_errors()
        cov = covar_results.extra_output

        niter = 10
        stats, accept, params = mcmc.get_draws(self.fit, cov, niter=niter)
        # try:
        #     assert (covar_results.parmaxes < params.std(1)).all()
        # except AssertionError:
        #     print 'covar: ', str(covar_results.parmaxes)
        #     print 'param: ', str(params.std(1))
        #     raise


if __name__ == '__main__':

    import sherpa.astro.sim as sim
    SherpaTest(sim).test(datadir="/data/scialg/testdata")