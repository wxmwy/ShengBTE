from pprint import pprint

import matplotlib.pyplot as plt

DGEMM_perf = [
    #482
    83.7
    #29413.07153,42096.60843,54697.74343,61806.01595,83220.31746,85850.39353,92395.31933,97943.4269,110592,102661.3847,106541.1363,119746.6299,138891.4244,138045.9475,137390.2556,135191.4057,144595.2729,137094.011,137389.2204,136677.4961,152489.3395,143789.9319,140967.6465,139580.277,147037.7008,138280.5745,134598.4802,148792.6667,168302.0485,163413.6975,160310.8686,159176.9958,166215.3447,156804.3442,153308.2967,158776.4398,171291.1657,163396.6638,160762.5583,158592.4484,164695.0477,157556.7492,155000.4397,168355.6507,188101.0016,180314.8479,177618.9692,175041.7989,181375.3974,172276.9743,167450.038,172187.8099,182781.8542,175843.5312,172682.6954,172335.3484,176595.3338,169665.3533,166502.2915,179785.6033,193078.8807,191897.1648,187149.5576,185332.913,190353.959,183270.6162,177714.84,182401.1268,192082.5972,185442.7825,180539.1005,180632.5668,185167.7212,180037.252,173841.1835,185632.384,200975.6483,197595.3801,193250.7653,191708.9178,197247.4608,190309.1206,185547.8426,188392.9874,198086.7207,191098.9576,187192.0108,186150.9344,190622.1946,185886.4206,181321.693,191476.8123,205699.1999,201701.676,198415.6809,196307.1357,201190.8252,195061.2752,191035.851,192690.0712,201182.5135,195237.3541,192500.1093,190344.4577,194692.4375,189282.085,186350.3048,195316.2939,208048.7346,205141.4114,201420.0771,199264.9003,203386.9302,198191.1279,194505.9677,196766.0051,203617.021,198245.1814,193886.7384,193623.1293,197243.151,192536.9953,188785.3197,197413.4853,209216.7396,206976.4175,203966.5698,201668.5177,204940.9434,199827.9863,196591.7161,197995.3253,204401.8898,199942.4437,196759.3388,195025.9813,197660.9045,193948.0425,190194.9434,199562.5658,210866.5914,208689.6895,205298.5297,203151.219,206326.5933,201700.6041,197220.8397,199887.1362,204911.677,201198.3241,197448.7943,196617.7203,196996.101,194381.03,190356.3409,201092.0658,211658.8567,210512.7338,205393.64,204831.9019,206564.7585,203060.2544,198088.5315,201559.7312,204168.0665,202399.2547,196969.0471,197880.3924,197184.6718,195502.2043,189521.9512,201484.2275,209759.3144,210463.1664,205147.206,205343.7147,205228.9018,202749.0656,195127.2699,201113.0837,200106.4966,201724.5499,196137.8347,197091.5329,195434.9676,195294.9103,188050.9655,201083.8905,211882.681,208364.9629,207432.5787,205558.8166,206673.0747,204569.516,205825.6241,203987.9267,207341.5785
]

DTRSV_perf = [
    482
    #570.155902,1177.011494,2276.679842,2807.883462,3124.237247,3572.43919,3489.776047,3720.254314,3791.205777,3873.945447,4275.2053,4408.119339,4469.421488,4871.929314,5009.457961,5439.236435,5510.604622,5638.134081,6044.508396,6015.390942,6127.243864,6548.837209,6580.529167,6926.637934,7054.829251,7009.65034,7244.57988,7298.393622,7252.259004,7679.104105,7723.600973,8104.872619,8076.715821,8106.447526,8441.908044,8413.396984,8402.447869,8807.080658,8722.043333,8997.995442,9054.034578,8955.779766,9268.215794,9255.676062,9215.795205,9541.857127,9506.7055,9796.032254,9788.801873,9767.040816,10132.94376,9980.161476,9829.666537,10140.09434,10063.93928,10313.40408,10427.13645,10376.55472,10636.93713,10591.12522,10527.15568,10826.05455,10778.69009,10971.43814,10953.16096,10893.61702,11142.00878,11081.43669,11037.97537,11333.29268,11256.82073,11441.34078,11435.11908,11369.58809,11602.30676,11579.90939,11506.38117,11691.87576,11556.71772,11766.54356,11739.32735,11753.5391,12026.54101,11980.50049,11911.30918,12165.3139,12073.87652,12278.74021,12287.55944,12186.92973,12402.05518,12346.64657,12252.21656,12496.20407,12457.08509,12687.06434,12627.14672,12560.93483,12816.29661,12717.6374
]

SPMV_AI = [
    1.42
    #0.0991883,0.0978661,0.0850436,0.098607,0.098607,0.0958967,0.0970855,0.0986889,0.0998498,0.0982512,0.0903348,0.098993,0.0964562,0.0995509,0.0984202,0.0940355,0.0954902,0.0979223,0.09204,0.0920917,0.0949505
]

SPMV_perf = [
    768.2
    #2.04381,1.28339,1.13244,1.84105,1.79029,0.283703,0.0293527,1.55788,2.19084,1.86182,1.03054,1.39991,1.60371,2.07569,1.79156,0.385395,0.306012,0.212026,0.158356,1.01184,0.110704
    #6.725116271,8.51560734,2.883717646,5.431340457,5.35112155,3.712256416,6.387275788,3.155108004,2.452960253,4.374496231,2.530690137,5.460768375,4.080438712,5.087949338,4.281257029,5.474005407,1.603910098
]

def frange(start, stop, step=1.0):
    f = start
    while f < stop:
        f += step
        yield f

Freq = 2.4 # GHz
Core = 14  # 32 cores
ThreadsPerCore = 1 # Estimated
Latency = 3 # 3 cycles per FADD instruction
SIMD_width = 2
SIMD_throughput = 1

#TLP = Freq * Core
#ILP = Freq * Core * max(1, ThreadsPerCore/Latency)
#SIMD= Freq * Core * SIMD_width / SIMD_throughput
#PEAK= 2 * SIMD # 2 SIMD parts in a core
PEAK = 1075

#pprint(TLP)
#pprint(ILP)
#pprint(SIMD)

max_flops = [PEAK,4700, 7000]
labels = ['CPU Theoretical Peak','P100 Theoretical Peak','V100 Theoretical Peak']
#max_flops = [PEAK,SIMD, ILP, TLP]
#labels = ['Theoretical Peak','SIMD', 'ILP', 'TLP']

# Plot configuration
height = 0.8
fig = plt.figure(frameon=False)
ax = fig.add_subplot(1, 1, 1)

yticks_labels = []
yticks = []
xticks_labels = []
xticks = [10.**i for i in range(-2, 6)]

ax.set_xlabel('Operational Intensity (Flops/Byte, log scale) ')
ax.set_ylabel('Performance (GFlops, log scale) ')
#ax.set_title('FTP Roofline Model')

# membw = [ 46.8114, 13.0187 ]
membw = [  76.8,732,900 ]

DGEMM_AI = [ 1.42 ]
#DGEMM_AI = [ x / 16 for x in [110592,124416,132096] ]
#DGEMM_AI = [ x / 16 for x in range(128, 120000, 32) ]
#DGEMM_perf = [ x / 1000 for x in DGEMM_perf ]
#SPMV_perf = [ x / 4 for x in SPMV_perf ]
DTRSV_AI = [ 1.42 ]
#DTRSV_perf = [ x / 1000 for x in DTRSV_perf ]

# Upper bound
x = list(frange(min(xticks), max(xticks), 0.01))
for i in range(0,1):
    ax.plot(x, [float(max_flops[0]) for x in x], label=labels[0], color='red',linestyle='--')
    bw = membw[0]
    ax.plot(x, [min(bw*x, float(max_flops[0])) for x in x],color='red', label=labels[0])
    
    ax.plot(x, [float(max_flops[1]) for x in x], label=labels[1], color='blue', linestyle='--')
    bw = membw[1]
    ax.plot(x, [min(bw*x, float(max_flops[1])) for x in x], color='blue',label=labels[1])

    ax.plot(x, [float(max_flops[2]) for x in x], label=labels[2], color='yellow', linestyle='--')
    bw = membw[2]
    ax.plot(x, [min(bw*x, float(max_flops[2])) for x in x], color='yellow',label=labels[2])

plt.text(150, max_flops[0]+10, labels[0], horizontalalignment='center')
plt.text(150, max_flops[1]+10, labels[1], horizontalalignment='center')
plt.text(150, max_flops[2]+10, labels[2], horizontalalignment='center')
#plt.text(0.1, max_flops[2]+5, "+" + labels[2], horizontalalignment='center')
#plt.text(0.1, max_flops[3]-25, "+" + labels[3], horizontalalignment='center')

#plt.text(54, 26, "HPL", horizontalalignment='center')
#plt.text(2, 0.16, "TRSV", horizontalalignment='center')
#plt.text(0.1, 0.06, "SpMV", horizontalalignment='center')

#plt.text(0.1, 30, "peak stream bandwidth", horizontalalignment='center', rotation=38)
#plt.text(0.4, 24, "+ Memory Affinity", horizontalalignment='center', rotation=40)

ylim_min, ylim_max = ax.get_ylim()

ax.plot(DGEMM_AI, DGEMM_perf, 'rx', markersize=4, markeredgewidth=1)
ax.plot(DTRSV_AI, DTRSV_perf, 'bx', markersize=4, markeredgewidth=1)
ax.plot(SPMV_AI, SPMV_perf, 'yx', markersize=4, markeredgewidth=1)

#ax.annotate("N Scales", xy=(3.5, 15), xytext=(3.5, 0.1),arrowprops=dict(arrowstyle="->"), horizontalalignment='center')

#arith_intensity = result['mem bottlenecks'][result['bottleneck level']]['arithmetic intensity']
#arith_intensity = [ 0.1, 0.2, 0.3, 0.4, 0.455729445, 0.5, 1, 2, 4 ]
#attain = []
#perf = []
#for i in arith_intensity:
#    perf.append(min(i*bw, max_flops[0]))
#ax.plot(arith_intensity, perf, 'r+', markersize=12, markeredgewidth=4)

# ax.tick_params(axis='y', which='both', left='off', right='off')
# ax.tick_params(axis='x', which='both', top='off')
#ax.set_xscale('log', basex=2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(min(xticks)+0.01, max(xticks))
# ax.set_yticks([perf, float(max_flops)])
#ax.set_xticks(xticks+arith_intensity)
#ax.grid(axis='x', alpha=0.7, linestyle='--')
fig.savefig('./pdf/roofline.pdf')
plt.show()
