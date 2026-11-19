import os
import sys, numpy as np, time
sys.path.insert(0, '.')
import dp
consts = {"pi": np.pi, "gamma": 0.5772156649015329, "delta": 4.669201609102991}
t=time.time(); dp.window_top(int(sys.argv[1]), sys.argv[2]=="1", consts); print("t",time.time()-t)
