import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

np.random.seed(1)

mus =  np.array([[0.2], [0.8]])
sigmas = np.array([[0.1], [0.1]]) ** 2
gmm = GaussianMixture(2)
gmm.means_ = mus
gmm.covars_ = sigmas
gmm.weights_ = np.array([0.5, 0.5])

#Fit the GMM with random data from the correspondent gaussians
gaus_samples_1 = np.random.normal(mus[0], sigmas[0], 10).reshape(10,1)
gaus_samples_2 = np.random.normal(mus[1], sigmas[1], 10).reshape(10,1)
fit_samples = np.concatenate((gaus_samples_1, gaus_samples_2))
print(fit_samples)
print(fit_samples.shape)
gmm.fit(fit_samples)

fig = plt.figure()
ax = fig.add_subplot(111)
x = np.linspace(0, 1, 1000).reshape(1000,1)
logprob = gmm.score_samples(x)
pdf = np.exp(logprob)
#print np.max(pdf) -> 19.8409464401 !?
ax.plot(x, pdf, '-k')
plt.show()


print(5**2)