import os

base_wiki_path = '/home/wkerzend/public_html/mcsnr/data/pages/moriagen.d'
base_image_path = '/home/wkerzend/public_html/mcsnr/data/media/moriagen.d'


class GeminiUtilDBMixin(object):

    geminiutil_session = None



class SNRGeminiTargetPlotting(object):
    def plot_cmd(cur_snr, candidates, figure=None):
        pass

def get_dict_uncertainty_string(stellar_params, param):
    param_value = stellar_params[param]
    param_uncertainy = stellar_params[param+'_uncertainty']

    if param_value is not None:
        param_value_str = '{0:.2f}'.format(param_value)
    else:
        param_value_str = 'none'

    if param_uncertainy is not None:
        param_uncertainy_str = '{0:.2f}'.format(param_uncertainy)
    else:
        param_uncertainy_str = 'none'

    return '{0} +/- {1}'.format(param_value_str, param_uncertainy_str)

class SNRGeminiTargetWiki(object):

    def make_wiki_page(self, image_width=600):
        snr_dir = os.path.join(base_wiki_path, '{0}.d'.format(self.snr.name))
        page_path = os.path.join(snr_dir, '{0}.txt'.format(self.snr.name))
        image_dir = os.path.join(base_image_path, '{0}.d'.format(self.snr.name))
        if not os.path.exists(snr_dir):
            os.system('mkdir -p {0}'.format(snr_dir))

        if not os.path.exists(image_dir):
            os.system('mkdir -p {0}'.format(image_dir))


        #os.system('mkdir -p %s' % wiki_image_path)
        """
        image_path = get_image_path(snr.snr.name)
        print "plot radio"
        plot_snr_simple(image_path, 'radioband.fits', cur_snr, figure=figure(1), save=os.path.join(wiki_image_path, 'simple_radioband.png'))
        print "plot rband"
        plot_snr_simple(image_path, 'rband.fits', cur_snr, figure=figure(1), save=os.path.join(wiki_image_path, 'simple_rband.png'))
        print "plot xray"
        plot_snr_simple(image_path, 'xrayband.fits', cur_snr, figure=figure(1), save=os.path.join(wiki_image_path, 'simple_xrayband.png'))
        print "plot color"
        plot_snr_color(image_path, cur_snr, save=os.path.join(wiki_image_path, 'simple_color.png'))

        candidates = session.query(Candidate).filter_by(snr_id=cur_snr.snr_id).all()
        if len(candidates) == 0:
            no_candidates = True
        else:
            no_candidates = False
        if not no_candidates:

            fitsfig = plot_snr_simple(image_path, 'rband.fits', cur_snr, figure=figure(1))
            overplot_candidates(fitsfig, candidates, save=os.path.join(wiki_image_path, 'simple_candidates.png'))
            fitsfig = plot_snr_simple(image_path, 'rband.fits', cur_snr, figure=figure(1) )
            overplot_candidates_label(fitsfig, candidates, save=os.path.join(wiki_image_path, 'simple_candidates_label.png'))

            fitsfig = plot_snr_simple(image_path, 'rband.fits', cur_snr, figure=figure(1))
            overplot_slits(fitsfig, cur_snr, candidates, save=os.path.join(wiki_image_path, 'simple_candidates_slit.png'))

            plot_cmd(cur_snr, candidates, figure=figure(1), save=os.path.join(wiki_image_path, 'simple_candidates_cmd.png'))

        """
        with file(page_path, 'w') as fh:
            fh.write('====== %s ======\n\n' % self.snr.name.upper())
            fh.write('**{0}**\n\n'.format(str(self).strip('<>')))

            fh.write('\n\n====== Candidates ======\n\n')
            fh.write('^ Name ^ MCPS ID ^ Coordinates ^ Teff ^ logg ^ [Fe/H] ^ vrad ^ vrot^\n')

            for candidate in self.candidates:
                columns = ['[[moriagen.d:{snrname}.d:star-{label}|'
                           'Star-{label}]]'.format(snrname=self.snr.name,
                                                   label=candidate.label)]

                columns.append('{0}'.format(candidate.mcps.id))
                columns.append('{0} {1}'.format(str(candidate.mcps.ra_sex),
                                                str(candidate.mcps.dec_sex)))

                stellar_params = candidate.get_combined_stellar_parameters()

                for param in ['teff', 'logg', 'feh', 'vrad', 'vrot']:
                        columns.append(get_dict_uncertainty_string(
                            stellar_params, param))

                fh.write('|{0}|\n'.format('|'.join(columns)))

            fh.write('\n\n====== Calibration Stars ======\n\n')

            fh.write('^ Name ^ MCPS ID ^ Coordinates ^ Teff ^ logg ^ [Fe/H] ^ vrad ^ vrot^\n')

            for candidate in self.calibration_stars:
                columns = ['[[moriagen.d:{snrname}.d:star-{label}|'
                           'Star-{label}]]'.format(snrname=self.snr.name,
                                                   label=candidate.label)]

                columns.append('{0}'.format(candidate.mcps.id))
                columns.append('{0} {1}'.format(str(candidate.mcps.ra_sex),
                                                str(candidate.mcps.dec_sex)))

                stellar_params = candidate.get_combined_stellar_parameters()

                for param in ['teff', 'logg', 'feh', 'vrad', 'vrot']:
                        columns.append(get_dict_uncertainty_string(
                            stellar_params, param))

                fh.write('|{0}|\n'.format('|'.join(columns)))



        os.system('chmod a+rw {0}'.format(page_path))
        os.system('chmod a+rwX {0}'.format(image_dir))

class CandidateWiki(object):
    def make_wiki_page(self, fig, isochrone_fitter):
        snr_dir = os.path.join(base_wiki_path, '{0}.d'.format(self.snr.name))
        page_path = os.path.join(snr_dir, 'star-{0}.txt'.format(self.label))
        image_dir = os.path.join(base_image_path, '{0}.d'.format(self.snr.name))

        with open(page_path, 'w') as fh:
            fh.write('==== Star {0} ====\n\n'.format(self.label))

            fh.write('<code>proj.session.query(ma.Candidate).filter_by(id={0}).one()</code>\n\n'.format(self.id))

            fh.write('^Spectrum^ Grating^ Teff ^ logg ^ [Fe/H] ^ vrad ^ vrot^\n')

            for i, mos_spectrum in enumerate(self.mos_spectra):
                grating_name = (mos_spectrum.slice.science_set.science.
                                instrument_setup.grating.name)

                stellar_params = self.preliminary_stellar_parameters[i]
                columns = ['Spectrum {0}'.format(i)]
                columns.append(grating_name)

                for param in ['teff', 'logg', 'feh', 'vrad', 'vrot']:
                        columns.append(getattr(stellar_params, param + '_str'))

                fh.write('|{0}|\n'.format('|'.join(columns)))

            fh.write('\n\n =====Plots====\n\n')
            fh.write('^ Fit ^ Isochrone^\n')

            #### Plotting fits and images ####
            for i, mos_spectrum in enumerate(self.mos_spectra):
                stellar_params = self.preliminary_stellar_parameters[i]
                fig.clf()
                fitimagename = 'star{0}-spectrum{1}-fit.png'.format(self.label, i)
                if len(mos_spectrum.wavelength) > 100 and \
                        (None not in [getattr(stellar_params, param)
                                      for param in ['teff', 'logg', 'feh',
                                                    'vrad', 'vrot']]):
                    ax = fig.add_subplot(111)
                    self.plot_fit(ax, i)


                else:
                    ax = fig.add_subplot(111)
                    ax.text(0.5, 0.5,'Bad extraction!!!',
                         horizontalalignment='center',
                         verticalalignment='center',
                         transform=ax.transAxes,
                         color='red',
                         fontsize=20)

                fig.savefig(os.path.join(image_dir, fitimagename))
                fig.clf()
                ax = fig.add_subplot(111)
                (age, metallicity), closest_isochrone = isochrone_fitter.find_closest_isochrone(
                    teff=stellar_params.teff,
                    logg=stellar_params.logg,
                    feh=stellar_params.feh)
                ax.plot(closest_isochrone['teff'], closest_isochrone['logg'])
                if (stellar_params.teff_uncertainty is not None) and \
                        (stellar_params.logg_uncertainty is not None):
                    ax.errorbar([stellar_params.teff],
                                [stellar_params.logg],
                                xerr=[stellar_params.teff_uncertainty],
                                yerr=[stellar_params.logg_uncertainty])

                ax.invert_xaxis()
                ax.invert_yaxis()

                isochroneimagename = 'star{0}-spectrum{1}-isochrone.png'.format(self.label, i)
                fig.savefig(os.path.join(image_dir, isochroneimagename))

                fh.write('|{{{{moriagen.d:{snrname}.d:{fitname}}}}}|'
                         '{{{{moriagen.d:{snrname}.d:{isochronename}}}}}|\n'.format(
                    snrname=self.snr.name, fitname=fitimagename, isochronename=isochroneimagename))
                fh.write('|Teff {0:.2f} logg {1:.2f} feh {2} age {3:.2f} z {4:.2f}||\n'.format(
                    stellar_params.teff, stellar_params.logg, stellar_params.feh, age, metallicity))

