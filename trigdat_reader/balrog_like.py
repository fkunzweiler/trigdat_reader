from threeML.plugins.DispersionSpectrumLike import DispersionSpectrumLike
from astromodels.functions.priors import Uniform_prior, Cosine_Prior


class BALROGLike(DispersionSpectrumLike):
    def __init__(self,
                 name,
                 observation,
                 background=None,
                 time=0,
                 free_position=True,
                 verbose=True):

        self._free_position = free_position

        #self._is_rsp_set = False

        super(BALROGLike, self).__init__(name, observation, background,
                                         verbose)

        # only on the start up
        self._rsp.set_time(time)

    def set_model(self, likelihoodModel):

        super(BALROGLike, self).set_model(likelihoodModel)
        if self._free_position:

            for key in self._like_model.point_sources.keys():

                self._like_model.point_sources[key].position.ra.free = True
                self._like_model.point_sources[key].position.dec.free = True

                self._like_model.point_sources[
                    key].position.ra.prior = Uniform_prior(
                        lower_bound=0., upper_bound=360)
                self._like_model.point_sources[
                    key].position.dec.prior = Cosine_Prior(
                        lower_bound=-90., upper_bound=90)

    def get_model(self):

        # Here we update the GBM drm parameters which creates and new DRM for that location
        # we should only be dealing with one source for GBM

        #if not self._is_rsp_set:

        for key in self._like_model.point_sources.keys():
            ra = self._like_model.point_sources[key].position.ra.value
            dec = self._like_model.point_sources[key].position.dec.value

        # update the location
        
        self._rsp.set_location(ra, dec)

        # if not self._free_position:

        #     # if we are not fitting for position, we only need to get the RA and DEC once

        #     self._is_rsp_set = True

        return super(BALROGLike, self).get_model()

    @classmethod
    def from_spectrumlike(cls, spectrum_like, time, free_position=True):

        return cls(spectrum_like.name, spectrum_like._observed_spectrum,
                   spectrum_like._background_spectrum, time, free_position,
                   spectrum_like._verbose)
