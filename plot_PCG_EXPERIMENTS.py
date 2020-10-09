import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mticker
import scipy.io
from matplotlib.lines import Line2D
def set_pars(mpl):
    mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath,bm,amsfonts}"]
    params = {'text.usetex': True,
              'font.family': 'serif',
              'font.size': 12,
              'legend.fontsize': 10,
              }
    mpl.rcParams.update(params)
    fig_par = {'dpi': 1000,
               'facecolor': 'w',
               'edgecolor': 'k',
               'figsize': (4, 3),
               'figsize3D': (4, 4),
               'pad_inches': 0.02,
               }

    return fig_par

n=49
x= np.arange(1,n+1)
diag= np.arange(1,n+1)
diagplus1=np.zeros(n)
diagminus1=np.zeros(n)
for i in np.arange(0,n):
    diagplus1[i] =(i+2)
    diagminus1[i] =(i)
    if (i)%7==0:
        diagminus1[i] = (i+7)*100

    if (i+1)%7==0:
         diagplus1[i] = (i-5)*100
downx=np.zeros(7)
downy=np.zeros(7)
for i in np.arange(7):
    downx[i] =i*7+1
    downy[i]=downx[i]+6

upx = np.zeros(7)
upy = np.zeros(7)
for i in np.arange(7):
    upx[i] = (i+ 1) * 7
    upy[i] = upx[i] - 6

cubex=np.array([0.5,7.5,7.5,0.5,0.5])
cubey=np.array([0.5,0.5,7.5,7.5,0.5])


fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\


plt.figure(num=None, figsize=parf['figsize3D'], dpi=parf['dpi'])

for i in np.arange(7):
    plt.plot(cubex+7*i,cubey+7*i,'--',color='blue',linewidth=0.5)

plt.plot(x, diag, 'ko', linewidth=0., marker='o', markeredgewidth=1,markersize=2,
         markerfacecolor='k')
plt.plot(x, diagplus1, 'ko', linewidth=0., marker='o', markeredgewidth=1,markersize=2,
         markerfacecolor='k')
plt.plot(x, diagminus1, 'ko', linewidth=0., marker='o', markeredgewidth=1,markersize=2,
         markerfacecolor='k')

plt.plot(downx, downy, 'rx', linewidth=0., marker='x', markeredgewidth=1,markersize=2,
         markerfacecolor='r')

plt.plot(upx, upy, 'rx', linewidth=0., marker='x', markeredgewidth=1,markersize=2,
         markerfacecolor='r')

# plt.xlabel(r'$\dagger = 10$')
plt.xlabel(r'(a)')
ax = plt.gca()
ax.set_xlim([0, 50])
ax.set_ylim([0,50])
plt.gca().invert_yaxis()
#plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))


colors = ['black', 'red', 'green']
lines = [Line2D([0], [0], linewidth=0., marker='o', markeredgewidth=1,markersize=3,
         markerfacecolor='k',color='k') ,
         Line2D([0], [0], linewidth=0., marker='x', markeredgewidth=1, markersize=3,
                markerfacecolor='r',color='r'),
         Line2D([0], [0], linestyle='--', linewidth=1., color='blue')
         ]
labels = ['non-zero el.', 'outlying non-zero el.', 'diagonal blocks']
plt.legend(lines, labels)

#plt.legend(loc='best')
fname = src + 'Sparsity{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot Sparsity')

plt.show()

breakpoint()


mat1 = scipy.io.loadmat('Spectra1.mat')
Ns = np.arange(mat1['SpectrumRed'].size)
S1 = mat1['SpectrumRed']
S2 = mat1['SpectrumBlue']
S3 = mat1['SpectrumBlack']

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.semilogy(Ns, S1, 'r|', linewidth=0., marker='|', markeredgewidth=0.5,markersize=2,
         markerfacecolor='None', label=r'$(\widehat{K}^{\textrm{ref}}_{0})^{-1}\widehat{K}$')
plt.semilogy(Ns, S2, 'bo', linewidth=0., marker='o', markeredgewidth=0.5,markersize=2,
         markerfacecolor='None', label=r'$(\widehat{K}^{\textrm{ref}}_{1})^{-1}\widehat{K}$')
plt.semilogy(Ns, S3, 'kx', linewidth=0., marker='x', markeredgewidth=0.5,markersize=2,
         markerfacecolor='None', label=r'$(\widehat{K}^{\textrm{ref}}_{2})^{-1}\widehat{K}$')

plt.xlabel(r'$\dagger = 1$')
plt.ylabel('eigenvalues')
ax = plt.gca()
ax.set_xlim([0, 444])
ax.set_ylim([1e-2,200])
#plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'Spectra1{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot Spectra1')

plt.show()



mat1 = scipy.io.loadmat('Spectra10.mat')
Ns = np.arange(mat1['SpectrumRed'].size)
S1 = mat1['SpectrumRed']
S2 = mat1['SpectrumBlue']
S3 = mat1['SpectrumBlack']

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.semilogy(Ns, S1, 'r|', linewidth=0., marker='|', markeredgewidth=0.5,markersize=2,
         markerfacecolor='None', label=r'$(\widehat{K}^{\textrm{ref}}_{0})^{-1}\widehat{K}$')
plt.semilogy(Ns, S2, 'bo', linewidth=0., marker='o', markeredgewidth=0.5,markersize=2,
         markerfacecolor='None', label=r'$(\widehat{K}^{\textrm{ref}}_{1})^{-1}\widehat{K}$')
plt.semilogy(Ns, S3, 'kx', linewidth=0., marker='x', markeredgewidth=0.5,markersize=2,
         markerfacecolor='None', label=r'$(\widehat{K}^{\textrm{ref}}_{2})^{-1}\widehat{K}$')

plt.xlabel(r'$\dagger = 10$')
plt.ylabel('eigenvalues')
ax = plt.gca()
ax.set_xlim([0, 444])
ax.set_ylim([1e-2,200])
#plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'Spectra10{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot Spectra10')

plt.show()




Ns = np.array(
    [1,	3,	5,	7,	9,	11,	13,	15,	17,	19])
S1 = np.array([21,	36,	45,	52,	58,	63,	68,	72,	76,	80])
S2 = np.array([21,	20,	20,	20,	20,	20,	19,	19,	19,	19])
S3 = np.array([23,	22,	22,	22,	22,	21,	21,	21,	21,	21])

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.plot(Ns, S1, 'r', linewidth=1., marker='|', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{0}$')
plt.plot(Ns, S2, 'b', linewidth=1., marker='o', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{1}$')
plt.plot(Ns, S3, 'k', linewidth=1., marker='x', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{2}$')

plt.xlabel(r' ration of eigenvalues - $\mu$ ')
plt.ylabel('Number of iteration \n to reach $10^{-6}$ residual norm')
ax = plt.gca()
ax.set_xlim([0, 20])
ax.set_ylim([0,100])
plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'stepsExp3{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot stepsExp3')

plt.show()










Ns = np.array(
    [9, 361, 9801, 16641, 26569, 83521, 114921, 154449, 203401, 263169, 335241,
     421201, 522729, 641601])
S1 = np.array([7, 16, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19])
S2 = np.array([7, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14])
S3 = np.array([6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, ])

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.plot(Ns, S1, 'r', linewidth=1., marker='|', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{0}$')
plt.plot(Ns, S2, 'b', linewidth=1., marker='o', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{1}$')
plt.plot(Ns, S3, 'k', linewidth=1., marker='x', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{2}$')
plt.xlabel(r' grid size - $ |\mathbb{Z}^d_{\bm{N}} |$ ')
plt.ylabel('Number of iteration \n to reach $10^{-6}$ residual norm')
# plt.title('Resid. norm evolution \n geom {}, contrast= {}')
ax = plt.gca()
ax.set_xlim([0, 0.7e6])
ax.set_ylim([0, 25])
plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'stepsExp1{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot stepsExp1')
plt.show()

Ns = np.array(
    [9, 361, 5329, 9801, 16641, 26569, 40401, 59049, 83521, 203401, 263169, 421201, 522729, 641601])
T1 = np.array([0.0164800000000000, 0.00687400000000000, 0.0460690000000000, 0.0446530000000000, 0.0636610000000000,
               0.214255000000000, 0.239456000000000, 0.260934000000000, 0.490800000000000, 1.15772700000000,
               1.51324700000000, 3.46832100000000, 4.28070500000000, 6.12351200000000])
T2 = np.array([0.00254800000000000, 0.00408000000000000, 0.0353160000000000, 0.0397200000000000, 0.0406570000000000,
               0.120192000000000, 0.159983000000000, 0.161327000000000, 0.369037000000000, 0.909900000000000,
               1.18773100000000, 2.28011300000000, 3.47365600000000, 4.49438000000000])
T3 = np.array([0.0194940000000000, 0.00273600000000000, 0.0210630000000000, 0.0213060000000000, 0.0262890000000000,
               0.0471810000000000, 0.0782300000000000, 0.118798000000000, 0.179674000000000, 0.493752000000000,
               0.646586000000000, 1.53128000000000, 1.97852700000000, 2.59077200000000])

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.plot(Ns, T1, 'r', linewidth=1., marker='|', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{0}$')
plt.plot(Ns, T2, 'b', linewidth=1., marker='o', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{1}$')
plt.plot(Ns, T3, 'k', linewidth=1., marker='x', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{2}$')
plt.xlabel(r' grid size - $ |\mathbb{Z}^d_{\bm{N}} |$ ')
plt.ylabel('computational time [s]')
# plt.title('Resid. norm evolution \n geom {}, contrast= {}')
ax = plt.gca()
ax.set_xlim([0, 0.7e6])
ax.set_ylim([0, 7])
plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'timeExp1{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot timeExp1')

plt.show()

Ns = np.array(
    [9,	2601,	40401,	154449,	335241,	522729,	641601])
S1 = np.array([3, 123, 128, 129, 130, 131, 130])
S2 = np.array([1, 79, 80, 79, 80, 80, 80])
S3 = np.array([1, 81, 83, 83, 84, 84, 84])

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.plot(Ns, S1, 'r', linewidth=1., marker='|', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{0}$')
plt.plot(Ns, S2, 'b', linewidth=1., marker='o', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{1}$')
plt.plot(Ns, S3, 'k', linewidth=1., marker='x', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{2}$')

plt.xlabel(r' grid size - $ |\mathbb{Z}^d_{\bm{N}} |$ ')
plt.ylabel('Number of iteration \n to reach $10^{-6}$ residual norm')
ax = plt.gca()
ax.set_xlim([0, 0.7e6])
ax.set_ylim([0,150])
plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'stepsExp2{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot stepsExp2')

plt.show()



Ns = np.array(
    [9,	2601,	40401,	154449,	335241,	522729,	641601])
T1 = np.array([0.0141670000000000,	0.123141000000000,	1.16926900000000,	4.97662100000000,	13.4557790000000,	22.5170350000000,	28.0394560000000])
T2 = np.array([0.000903000000000000,	0.0735950000000000,	0.705271000000000,	3.04815200000000,	8.55199800000000,	13.8671200000000,	17.2844670000000])
T3 = np.array([0.00758200000000000,	0.0727510000000000,	0.756571000000000,	3.28678300000000,	9.06458800000000,	14.3490620000000,	17.5798820000000])

fig = plt.figure()
parf = set_pars(mpl)
src = './experiments/figures/'  # source folder\

plt.figure(num=None, figsize=parf['figsize'], dpi=parf['dpi'])

plt.plot(Ns, T1, 'r', linewidth=1., marker='|', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{0}$')
plt.plot(Ns, T2, 'b', linewidth=1., marker='o', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{1}$')
plt.plot(Ns, T3, 'k', linewidth=1., marker='x', markeredgewidth=1,
         markerfacecolor='None', label=r'$\widehat{K}^{\textrm{ref}}_{2}$')
plt.xlabel(r' grid size - $ |\mathbb{Z}^d_{\bm{N}} |$ ')
plt.ylabel('computational time [s]')
# plt.title('Resid. norm evolution \n geom {}, contrast= {}')
ax = plt.gca()
ax.set_xlim([0, 0.7e6])
ax.set_ylim([0, 30])
plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))

plt.legend(loc='best')
fname = src + 'timeExp2{}'.format('.pdf')
print(('create figure: {}'.format(fname)))
plt.savefig(fname, dpi=parf['dpi'], pad_inches=parf['pad_inches'], bbox_inches='tight')
print('END plot timeExp2')

plt.show()



