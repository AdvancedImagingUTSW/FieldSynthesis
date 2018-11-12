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

# Create F, the illumination pattern
F_hat = createAnnulus()
F_hat = ft.ifftshift(F_hat)
F = ft.ifft2(F_hat)
F = ft.fftshift(F)
# This is the illumination intensity pattern
Fsqmod = np.real(F*np.conj(F))

plt.figure()
plt.title('F')
plt.imshow(Fsqmod, cmap='plasma')
plt.show(block=False)

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

plt.figure()
plt.title('L')
plt.imshow(L, cmap='plasma')
plt.show(block=False)

plt.figure()
plt.title('L_hat')
plt.imshow(np.abs(L_hat), cmap='plasma')
plt.show(block=False)

# Manually scan by shifting Fsqmod and multiplying by Lsqmod
scanned = np.zeros(Fsqmod.shape)

for x in range(np.size(Fsqmod,1)):
    scanned = scanned + np.roll(Fsqmod,x-center,1)*Lsqmod[center,x]

plt.figure()
plt.title('Scanned')
plt.imshow(scanned, cmap='plasma')
plt.show(block=False)

# Manually scanning is a convolution operation
# There are potentially boundary effects here
convolved = sig.fftconvolve(Fsqmod,Lsqmod,'same')

plt.figure()
plt.title('Convolved')
plt.imshow(convolved, cmap='plasma')
plt.show(block=False)

# This manual implementation of Fourier transform based convolution
# actually does circular convolution
convolvedft = ft.fftshift(ft.fft2(ft.ifft2(ft.ifftshift(Fsqmod)) *ft.ifft2(ft.ifftshift(Lsqmod))))
convolvedft = np.real(convolvedft)

plt.figure()
plt.title('Convolved FT')
plt.imshow(convolvedft, cmap='plasma')
plt.show(block=False)

# Do the Field Synthesis method of performing a line scan at the back focal plane
fieldSynthesis = np.zeros_like(Fsqmod)

for a in range(fieldSynthesis.shape[1]):
    # Instaneous scan in frequency space
    temp_hat = F_hat * np.roll(L_hat,a-fieldSynthesis.shape[1]//2,1)
    # Instaneous scan in object space
    temp = ft.fftshift( ft.fft2( ft.ifftshift(temp_hat) ) )
    # Incoherent summing of the intensities
    fieldSynthesis = fieldSynthesis + np.abs(temp)**2

plt.figure()
plt.title('Field Synthesis')
plt.imshow(fieldSynthesis, cmap='plasma')
plt.show(block=True)
