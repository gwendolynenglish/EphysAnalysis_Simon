import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = [
[2.917,
2.916,
2.919,
2.917,
2.916],

[2.919,
2.917,
2.916,
2.917,
2.916,
2.917,
2.916],

[2.908,
2.908,
2.908,
2.909,
2.909,
2.911,],

[2.852,
2.967,
2.850,
2.963,
2.941,],

[2.761,
2.786,
2.775,
2.764,
2.751,
2.770,],

[2.708,
2.678,
2.694,
2.627,
2.721,
2.674,
2.697,
2.690,
2.721,
2.718,],

[
2.534,
2.526,
2.521,
2.531,
2.579,
2.559,
2.571,
2.515,],

[2.278,
2.277,
2.234,
2.226,
2.290,
2.270,
2.207,],

[2.108,
2.104,
2.114,
2.124,
2.108,
2.135,
2.127,
2.117,],

[1.937,
1.902,
1.911,
1.929,
1.924,
1.935,
1.911,
1.924,
1.914,
1.913,],

[1.777,
1.783,
1.772,
1.827,
1.815,
1.791,
1.827,
1.822,
],
]




new_dat = [
[2.920,
2.863,
2.971,
2.967,
2.863,
2.968,
2.896,],

[2.904,
2.834,
2.930,
2.852,
2.837,
2.936,
2.895,
2.903,
],


[2.840,
2.872,
2.858,
2.809,
2.879,
2.896,
2.876,],


[2.836,
2.834,
2.844,
2.850,
2.836,
2.818,
2.756],

[2.735,
2.789,
2.761,
2.786,
2.751,
2.783,
2.778,
2.764],


[2.689,
2.716,
2.678,
2.676,
2.687,
2.724,
2.724,
2.674],]

print(new_dat)
new_avgs = [np.mean(arr) for arr in new_dat]

new_xdat = np.arange(2000, 7001, 1000)

plt.loglog(new_xdat, new_avgs)

plt.grid(which='minor')


avgs = [np.mean(arr) for arr in data]
print(avgs)


x_dat = [100, 200, 500, 2000, 5000, 7000, 9000, 13000, 15000, 18000, 20000]







plt.loglog(x_dat, avgs)
ax = plt.gca()

ax.set_title('RC low pass filter')


ax.set_xlabel('frequency in Hz')
ax.set_ylabel('attentuation ratio in dB')



plt.show()