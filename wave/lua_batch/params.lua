-- Example script: time advecting Gaussian bump on 1e6 points

sim = wave1d.new{ n=1000000, c=1.0 }

t1 = wtime()
sim:set_frame(function(x) return math.exp(-25*(x-0.5-sim.c*sim.dt)^2) end, -1)
sim:set_frame(function(x) return math.exp(-25*(x-0.5)^2) end, 0)
t2 = wtime()
sim:step(1000)
t3 = wtime()

flops = 5.0 * sim.tidx * sim.n
print("Initialization: ", t2-t1)
print("Stepping: ", t3-t2)
print("Effective Gflop/s:", flops/(t3-t2)/1e9)
