import numpy as np
import matplotlib.pyplot as plt
import pyfits as p
import kplr
from scipy.ndimage.measurements import label
import mahotas as mh
from scipy.ndimage.measurements import center_of_mass
from scipy.interpolate import splrep, splev
import os
import files
import matplotlib.image as mpimg

fmt = ['ko', 'rD', 'b^', 'gs']




class star(object):
    """
    Routines to perform photometry on the Kepler full frame images
    
    Ben Montet, 2017
    https://github.com/benmontet/f3
    See repository for a tutorial
    """
    
    
    def __init__(self, kic, ffi_dir=None):
        self.obs_filenames = files.ffilist
        self.kic = kic
        if ffi_dir == None:
            self.ffi_dir = 'ffidata/'
        else:
            self.ffi_dir = ffi_dir
        
        dir = os.path.dirname(__file__)
        try:
            obs_info = os.path.join(dir, '../obs_info.txt')
            self.times, self.qs, self.year = np.loadtxt(obs_info, unpack=True)
        except:
            obs_info = os.path.join(dir, 'obs_info.txt')
            self.times, self.qs, self.year = np.loadtxt(obs_info, unpack=True)


 
         
    def make_postcard(self, npix=300, shape=(1070, 1132), buffer_size=15):
        """
        Develop a "postcard" region, npix by npix in size, around the target star.
        Other stars in this postcard will be used as possible reference stars.
        """
        
        source = self.kic
        client = kplr.API()
        targ = client.target(source)
        
        channel = [targ.params['Channel_0'], targ.params['Channel_1'], 
                    targ.params['Channel_2'], targ.params['Channel_3']]

        col = [targ.params['Column_0'], targ.params['Column_1'], 
                    targ.params['Column_2'], targ.params['Column_3']] 

        row = [targ.params['Row_0'], targ.params['Row_1'], 
                    targ.params['Row_2'], targ.params['Row_3']] 
        

        if None in row:
            raise ValueError('Star not on detector all quarters!')

        if None in col:
            raise ValueError('Star not on detector all quarters!')
            
        center = np.array([npix/2, npix/2])
        
        
        # If star close to edge, shift frame so that we have the full npix by npix
        # In this case, postcard will not be centered on target star
        if (np.min(col) < npix/2):
            jump = npix/2 - np.min(col) + buffer_size
            col += jump
            center[1] -= jump
            
        if (np.min(row) < npix/2):
            jump = npix/2 - np.min(row) + buffer_size
            row += jump
            center[0] -= jump
            
        if (np.max(row) > shape[0] - npix/2):
            jump = shape[0]-npix/2 - np.max(row) - buffer_size
            row += jump
            center[0] -= jump
            
        if (np.max(col) > shape[1] - npix/2):
            jump = shape[1]-npix/2 - np.max(col) - buffer_size
            col += jump
            center[1] -= jump
        
        fin_arr = np.zeros((len(self.times), npix, npix))

        for icount, iname in enumerate(self.obs_filenames): 
            a = p.open(self.ffi_dir+iname)

                        
            quarter = a[0].header['quarter']

            if int(quarter) == 0:
                season = 3
            else:
                season = (int(quarter) - 2) % 4

            #season_arr[icount] = season
            img = a[channel[season]].data
            img -= np.median(img)
            
            
            ymin = int(max([int(row[season])-npix/2,0]))
            ymax = int(min([int(row[season])+npix/2,img.shape[0]]))
            xmin = int(max([int(col[season])-npix/2,0]))
            xmax = int(min([int(col[season])+npix/2,img.shape[1]]))

            pimg = img[ymin:ymax,xmin:xmax]
            fin_arr[icount,:,:] = pimg
        
        self.postcard = fin_arr
        self.integrated_postcard = np.sum(self.postcard, axis=0)
        self.center = center
        
    
    def onclick(self, event):
        global ix, iy
        if event.xdata is not None and event.ydata is not None and event.xdata < 30:
            ix, iy = int(round(event.xdata)), int(round(event.ydata))
            global coords 
            self.coordsx.append(ix)
            self.coordsy.append(iy)
            
    
    def mini_lc(self, s1, s2, g, wh, factr):
        data_new = np.roll(self.postcard, s1, axis=1)
        data_new = np.roll(data_new, s2, axis=2)

        numer_pix = data_new[:,self.targets == 1]
        numer = np.sum(numer_pix, axis=1)

        factr[g] = numer[g] / self.reference_flux[g]
        factr[g] /= np.median(factr[g[wh]])
        return factr[g]


    def do_rolltest(self, g, wh):
        """
        Test for shifts in particular orientations by differential velocity aberrations.
        Allow for the possibility that the aperture could be shifted by up to one pixel
        as judged by the observed scatter
        """
        stdval_b = 1.0
        best = np.zeros(2)
        
        factr = np.zeros_like(self.reference_flux)
        
        for s1 in np.array([-1, 0, 1]):
            for s2 in np.array([-1, 0, 1]):
                factr[g] = self.mini_lc(s1, s2, g, wh, factr)

                stdval = np.zeros(4)
                for ij in xrange(4):
                    thisyear = np.where(self.year[g] == ij)[0]
                    if len(thisyear) >= 3:
                        stdval[ij] = np.std(factr[g][thisyear]/np.median(factr[g][thisyear]))

                stest = np.max(stdval[stdval != 0.0])
                fitline = np.polyfit(self.times[g][wh], factr[g][wh], 1)
                stest = np.max([np.std(factr[g][wh]/(fitline[0]*self.times[g][wh]+fitline[1])), 0.0001])

                if stest < stdval_b:
                    best = [s1, s2]

                    stdval_b = stest + 0.0
        return best

        
    
    def find_other_sources(self, plot_flag = False, edge_lim = 0.015, min_val = 5000, 
                          ntargets = 250, extend_region_size=3, remove_excess=4,
                          plot_window=15):
        """
        Identify apertures for all sources on the postcard, both for the 
        target and potential reference stars
        """
        j,i = self.center
        
        region = self.integrated_postcard + 0.0
        if plot_flag == True:
            ff = plt.imshow(self.integrated_postcard, interpolation='nearest', cmap='gray', 
                            vmax = np.percentile(region, 99.6))
            plt.colorbar(ff)
            plt.show()
        targets = np.zeros_like(self.integrated_postcard)
        sizeimg = np.shape(targets)[0]
        
        jj = j + 0
        ii = i + 0
        
        edge = edge_lim
        
        lim = max(min_val, self.integrated_postcard[j,i]*edge)
        
        maxpt = np.percentile(self.integrated_postcard, 94)
        
        bin_img = (region > lim)
        lab_img, n_features = label(bin_img)
        
        key_targ = (lab_img == (lab_img[j,i]))
        tot = np.sum(key_targ)
    
        targets[key_targ] = 1
        region[key_targ] = 0.0
        

        
        lim = np.zeros(ntargets)
        for peaks in xrange(1,ntargets):
            k = np.argmax(region)
            j,i = np.unravel_index(k, region.shape)
            lim[peaks] = max(maxpt, edge*region[j,i])


            bin_img = (region >= lim[peaks])
            lab_img, n_features = label(bin_img)

            key_targ = (lab_img == (lab_img[j,i]))
            targets[key_targ] = peaks + 1
            region[key_targ] = 0.0
        
        lab_img, n_features = label(targets)
        for i in xrange(1, ntargets+1):
            for j in xrange(extend_region_size):
                border= mh.labeled.border(targets, 0, i)

                targets[border*(region < (10)*lim[peaks])] = i
        for i in xrange(2, ntargets+1):
            for j in xrange(2, ntargets+1):
                if i != j:
                    border = mh.labeled.border(targets, i, j)
                    if np.sum(border) != 0:
                        targets[targets == j] = i
    
        targets = mh.labeled.remove_bordering(targets)
        for k in xrange(remove_excess):
            for i in xrange(ntargets):
                if np.sum(self.integrated_postcard[targets == i]) < 0.01:
                    targets[targets > i] -= 1

        self.targets = targets
        
        if plot_flag == True:
            plt.imshow(self.targets, interpolation='nearest')
            plt.show()

            plt.imshow(((targets == 1)*self.integrated_postcard + (targets == 1)*100000)
                       [jj-plot_window:jj+plot_window,ii-plot_window:ii+plot_window], 
                       interpolation='nearest', cmap='gray', vmax=np.percentile(self.integrated_postcard, 99.6))
            plt.show()

            plt.imshow((np.ceil(targets/100.0)*self.integrated_postcard+np.ceil(targets/500.0)*3500000), 
                       interpolation='nearest', cmap='gray', vmax=np.percentile(self.integrated_postcard, 99.99))
            plt.show()
            
    def do_photometry(self):
        """
        Actually calculate photometry
        """
        
        std_f = np.zeros(4)
        data_save = np.zeros_like(self.postcard)
        self.obs_flux = np.zeros_like(self.reference_flux)


        for i in xrange(4):
            g = np.where(self.qs == i)[0]
            wh = np.where(self.times[g] > 54947)

            data_save[g] = np.roll(self.postcard[g], int(self.roll_best[i,0]), axis=1)
            data_save[g] = np.roll(data_save[g], int(self.roll_best[i,1]), axis=2)

            self.target_flux_pixels = data_save[:,self.targets == 1]
            self.target_flux = np.sum(self.target_flux_pixels, axis=1)
            
            self.obs_flux[g] = self.target_flux[g] / self.reference_flux[g]
            self.obs_flux[g] /= np.median(self.obs_flux[g[wh]])
            
            fitline = np.polyfit(self.times[g][wh], self.obs_flux[g][wh], 1)
            std_f[i] = np.max([np.std(self.obs_flux[g][wh]/(fitline[0]*self.times[g][wh]+fitline[1])), 0.001])
        
        self.flux_uncert = std_f
        
    def generate_panel(self, img):
        """
        Create the figure viewed in adjust_aperture
        """
        plt.figure(figsize=(14,6))
        ax = plt.gca()
        fig = plt.gcf()
        plt.subplot(122)
        
        
        data_save = np.zeros_like(self.postcard)
        self.roll_best = np.zeros((4,2))
        
        for i in xrange(4):
            g = np.where(self.qs == i)[0]
            wh = np.where(self.times[g] > 54947)

            self.roll_best[i] = self.do_rolltest(g, wh)
            
        self.do_photometry()
        for i in xrange(4):
            g = np.where(self.qs == i)[0]
            plt.errorbar(self.times[g], self.obs_flux[g], yerr=self.flux_uncert[i], fmt=fmt[i])
            
        plt.xlabel('Time', fontsize=20)
        plt.ylabel('Relative Flux', fontsize=20)
 
        
        plt.subplot(121)
        implot = plt.imshow(img, interpolation='nearest', cmap='gray', vmin=98000*52, vmax=104000*52)
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        
        plt.show(block=True)

            
        
    def adjust_aperture(self, image_region=15, ignore_bright=0):
        """
        Develop a panel showing the current aperture and the light curve as judged from that aperture.
        Clicking on individual pixels on the aperture will toggle those pixels on or off into the
        aperture (which will be updated after closing the plot).
        Clicking on the 0th row or column will turn off all pixels in that column or row, respectively.
        Will iterate continuously until the figure is closed without updating any pixels.
        """
        self.ignore_bright = ignore_bright
        self.calc_fluxes()
        
        self.coordsx = []
        self.coordsy = []
        
        jj, ii = self.center
        plt.ion()
        img = np.sum(((self.targets == 1)*self.postcard + (self.targets == 1)*100000)
                     [:,jj-image_region:jj+image_region,ii-image_region:ii+image_region], axis=0)
        self.generate_panel(img)
        
        while len(self.coordsx) != 0:
            for i in xrange(len(self.coordsx)):
                if self.targets[self.coordsy[i]+jj-image_region,self.coordsx[i]+ii-image_region] != 1:
                    self.targets[self.coordsy[i]+jj-image_region,self.coordsx[i]+ii-image_region] = 1
                elif self.targets[self.coordsy[i]+jj-image_region,self.coordsx[i]+ii-image_region] == 1:
                    self.targets[self.coordsy[i]+jj-image_region,self.coordsx[i]+ii-image_region] = 0
                if self.coordsy[i] == 0:
                    thiscol = np.where(self.targets[:,self.coordsx[i]+ii-image_region] == 1)
                    self.targets[thiscol,self.coordsx[i]+ii-image_region] = 0
                if self.coordsx[i] == 0:
                    thiscol = np.where(self.targets[self.coordsy[i]+jj-image_region,:] == 1)
                    self.targets[self.coordsy[i]+jj-image_region, thiscol] = 0


            self.coordsx = []
            self.coordsy = []
            img = np.sum(((self.targets == 1)*self.postcard + 
                          (self.targets == 1)*100000)[:,jj-image_region:jj+image_region,ii-image_region:ii+image_region],
                         axis=0)
            self.generate_panel(img)

 
    def data_for_target(self, do_roll=True):
        """
        Determine the normalized photometry, accounting for effects shared by reference stars
        """
        self.calc_fluxes()
        self.roll_best = np.zeros((4,2))
        
        if do_roll == True:
            for i in xrange(4):
                g = np.where(self.qs == i)[0]
                wh = np.where(self.times[g] > 54947)

                self.roll_best[i] = self.do_rolltest(g, wh)
                
        self.do_photometry()

    
    def calc_fluxes(self, min_flux = 5000, outlier_iterations=5,
                       max_outlier_obs=4, outlier_limit=1.7):
        """
        Determine the suitable reference stars, and then the total flux in those stars and 
        in the target star in each epoch
        """
        
        jj, ii = self.center
        

        
        numer = np.zeros(len(self.times))
        denom = np.zeros(len(self.times))
        factr = np.zeros(len(self.times))
        
        numer_pix = self.postcard[:,self.targets == 1]
        numer = np.sum(numer_pix, axis=1)
        
        tar_vals = np.zeros((len(self.times), int(np.max(self.targets)+1-2-self.ignore_bright)))
        
        for i in xrange(2+self.ignore_bright,int(np.max(self.targets)+1)):
            tval = np.sum(self.postcard[:,self.targets == i], axis=1)
            #denom += tval/np.median(tval)
            tar_vals[:,i-2-self.ignore_bright] = tval #/ np.median(tval)
            
        for i in xrange(len(self.obs_filenames)):
            if np.max(tar_vals[i]) < 5000:
                tar_vals[self.qs == self.qs[i]] = 0.0

        all_tar = np.zeros((len(self.times), int(np.max(self.targets)-self.ignore_bright)))
        all_tar[:,0] = numer
        all_tar[:,1:] = tar_vals
        
        self.photometry_array = all_tar
        
        for i in xrange(len(tar_vals[0])):
            for j in xrange(4):
                g = np.where(self.qs == j)[0]  
                tar_vals[g,i] /= (np.median(tar_vals[g,i])+1e-15)
                

    
        for i in xrange(outlier_iterations):
            nonzeros = np.where(tar_vals[0,:] != 0)[0]
            med = np.median(tar_vals[:,nonzeros], axis=1)
            std = np.std(tar_vals[:,nonzeros], axis=1)

            if np.sum(tar_vals) != 0.0:
                tar_vals_old = tar_vals + 0.0

            for k in xrange(len(tar_vals[0])):
                h = np.where((np.abs(med-tar_vals[:,k])/std) > outlier_limit)[0]
                if len(h) >= max_outlier_obs:
                    tar_vals[:,k] = 0

        if np.sum(tar_vals) == 0.0:
            tar_vals = tar_vals_old + 0.0

        denom = np.sum(tar_vals, axis=1)
        self.target_flux_pixels = numer_pix
        self.reference_flux = denom
        
    def calc_centroids(self):
        """
        Identify the centroid positions for the target star at all epochs
        """
        self.cm = np.zeros((len(self.postcard), 2))
        for i in xrange(len(self.postcard)):
            target = self.postcard[i]
            target[self.targets != 1] = 0.0
            self.cm[i] = center_of_mass(target)
            
    def define_spotsignal(self):
        """
        Identify the "expected" flux value at the time of each observation based on the 
        Kepler long-cadence data.
        """
        client = kplr.API()
        star = client.star(self.kic)

        lcs = star.get_light_curves(short_cadence=False)
        time, flux, ferr, qual = [], [], [], []
        for lc in lcs:
            with lc.open() as f:
                hdu_data = f[1].data
                time.append(hdu_data["time"])
                flux.append(hdu_data["pdcsap_flux"])
                ferr.append(hdu_data["pdcsap_flux_err"])
                qual.append(hdu_data["sap_quality"])
            tout = np.array([])
            fout = np.array([])
            eout = np.array([])
            for i in xrange(len(flux)):
                t = time[i][qual[i] == 0]
                f = flux[i][qual[i] == 0]
                e = ferr[i][qual[i] == 0]

                t = t[np.isfinite(f)]
                e = e[np.isfinite(f)]
                f = f[np.isfinite(f)]

                e /= np.median(f)
                f /= np.median(f)
                tout = np.append(tout, t[50:]+54833)
                fout = np.append(fout, f[50:])
                eout = np.append(eout, e[50:])

            self.spot_signal = np.zeros(52)

            for i in xrange(len(self.times)):
                if self.times[i] < 55000:
                    self.spot_signal[i] = 1.0
                else:
                    self.spot_signal[i] = fout[np.abs(self.times[i] - tout) == np.min(np.abs(self.times[i] - tout))]

