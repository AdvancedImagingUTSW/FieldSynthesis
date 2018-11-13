import numpy as np
import scipy as sci
import scipy.fftpack as ft
import scipy.signal as sig
from scipy.stats import norm
import matplotlib.pyplot as plt

''' Field Synthesis
Python-based demonstration

Mark Kittisopikul
Goldman Lab
Northwestern University

November 2018
'''

def createAnnulus(n=256, r=32, w=4):
    ''' createAnnulus - create a ring-like structure
    INPUT
    n - size of square array or vector
    r - radius of the ring
    w - width of the ring
    OUTPUT
    an array n x n
    '''
    if np.isscalar(n):
        v = np.arange(n)
        v = v - np.floor(n/2)
    else:
        v = n

    y,x = np.meshgrid(v,v)
    q = np.hypot(x,y)
    annulus = abs(q-r) < w
    return annulus

plt.rc('text', usetex=True)
fig, ax = plt.subplots(2,4,sharey=True,sharex=True,figsize=(16,9))

# Create F, the illumination pattern
F_hat = createAnnulus()
F_hat = ft.ifftshift(F_hat)
F = ft.ifft2(F_hat)
F = ft.fftshift(F)
# This is the illumination intensity pattern
Fsqmod = np.real(F*np.conj(F))

#plt.figure()
#plt.title('F')
#plt.imshow(Fsqmod, cmap='plasma')
#plt.show(block=False)
ax[0,0].imshow(Fsqmod, cmap='plasma')
ax[0,0].set_title('F(x,z)')

# Create L, the scan profile
L = np.zeros_like(Fsqmod)
center = L.shape[1]//2
sigma = 30
L[center,:] = norm.pdf(np.arange(-center,center),0,sigma)
# L[L.shape[1]//2,:] = 1
# The square modulus of L is the object space
Lsqmod = L*np.conj(L)
# This is the line scan profile used in Field Synthesis
L_hat = ft.fftshift(ft.fft2(ft.ifftshift(L)))

ax[0,1].imshow(L, cmap='plasma')
ax[0,1].set_title('$ L(x)\delta(z) $')

ax[0,2].imshow(Lsqmod, cmap='plasma')
ax[0,2].set_title('$ |L(x)\delta(z)|^2 $')

ax[0,3].imshow(np.abs(L_hat), cmap='plasma')
ax[0,3].set_title('$\hat{L}(k_x) $')

# Manually scan by shifting Fsqmod and multiplying by Lsqmod
scanned = np.zeros(Fsqmod.shape)

for x in range(np.size(Fsqmod,1)):
    scanned = scanned + np.roll(Fsqmod,x-center,1)*Lsqmod[center,x]

ax[1,0].imshow(scanned, cmap='plasma')
ax[1,0].set_title('Scanned: $ \sum_{x\'} |F(x\',z)|^2|L(x-x\',z)|^2 $')

# Manually scanning is a convolution operation
# There are potentially boundary effects here
convolved = sig.fftconvolve(Fsqmod,Lsqmod,'same')

ax[1,1].imshow(convolved, cmap='plasma')
ax[1,1].set_title('Convolved: $ |F(x,z)|^2 ** |L(x,z)\delta(z)|^2 $')

# This manual implementation of Fourier transform based convolution
# actually does circular convolution
convolvedft = ft.fftshift(ft.fft2(ft.ifft2(ft.ifftshift(Fsqmod)) *ft.ifft2(ft.ifftshift(Lsqmod))))
convolvedft = np.real(convolvedft)

ax[1,2].imshow(convolvedft, cmap='plasma')
ax[1,2].set_title(r'Convolved FT: $ \mathcal{F}^{-1} \{ \mathcal{F}\{|F|^2\} \mathcal{F}\{|L|^2\} \} $')

# Do the Field Synthesis method of performing a line scan at the back focal plane
fieldSynthesis = np.zeros_like(Fsqmod)

for a in range(fieldSynthesis.shape[1]):
    # Instaneous scan in frequency space
    temp_hat = F_hat * np.roll(L_hat,a-fieldSynthesis.shape[1]//2,1)
    # Instaneous scan in object space
    temp = ft.fftshift( ft.fft2( ft.ifftshift(temp_hat) ) )
    # Incoherent summing of the intensities
    fieldSynthesis = fieldSynthesis + np.abs(temp)**2

ax[1,3].imshow(fieldSynthesis, cmap='plasma')
ax[1,3].set_title('Field Synthesis: $ \sum_a |\hat{F}(k_x,k_z)\hat{L}(k_x-a)|^2 $')

plt.show(block=True)
