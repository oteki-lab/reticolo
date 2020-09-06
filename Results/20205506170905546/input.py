""" reticolo input parameters """
asymmetry = True   # True: Combination only when the x-axis and the y-axis have the same value

npoints = 21        # point number of wavelength
lambdamin = 0.4     # min wavelength
lambdamax = 1.2     # max wavelength

## x
# Number of Fourier terms
Mx = 0
# period of nano structure
diam_x = 0.215
# height of nano structure
height_nanostructure = 0.4
# period of back grating
period_x = [2.4, 2.6]
# height of back grating
height_backgrating = 0.9

My = Mx                 # Number of Fourier terms in y
period_y = 4.8          # period in y
diam_y = diam_x
