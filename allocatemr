
# This file "allocatemr" contains the function: allocatemr().
# Version: 2017-10-28
# Author: Jaap de Gruijter

# CONTENT:
# 1. Systematic sampling from the grid;
# 2. Application of Ospats to the sample;
# 3. Allocation to sample-strata;
# 4. Sample sizes and Relative Standard Error;
# 5. Stratified random sampling
# 6. Output of final results.

# For process monitoring uncomment println(...) lines

function allocatemr()
println("-------------------- START FUNCTION ALLOCATEMR -----------------------")

##### SECTION 1. SYSTEMATIC SAMPLING FROM THE GRID WITH RANDOM STARTING POINT
inrow = 1:in
start = rand(rng, inrow)[1]
sample = start:in:N
n = length(sample)
println("sample size = ", n)
x_s = x[sample]
y_s = y[sample]
z_pred_s = z_pred[sample]
s2_s = s2[sample]

##### SECTION 2. APPLICATION OF OSPATS METHOD TO THE SAMPLE
##### Subsection 2.1 Calculate n x n matrix of generalized distances
# println("----- Calculating matrix of generalized distances --")
Lags = zeros(n,n)
Z2 = zeros(n,n)
SS = zeros(n,n)
for i = 1:(n-1)
for j = (i+1):n
    Lags[i,j] = sqrt((x_s[i] - x_s[j])^2 + (y_s[i] - y_s[j])^2)
    Z2[i,j] = (z_pred_s[i] - z_pred_s[j])^2
    SS[i,j] = s2_s[i] + s2_s[j]
end
end
Lags = Lags + Lags'
Z2 = Z2 + Z2'
SS = SS + SS'
Z2 = Z2/R2
  # Cov = 0.5*SS.*exp(-3*Lags/range)   explanation: 0.5*SS = Covmax
  # d2 = Z2 + SS - 2*Cov
d2 = Z2 + SS - SS.*exp(-3*Lags/range)

TOTd2 = sum(d2)/2
ObarH1 = sqrt(TOTd2)/n
# println("ObarH1 = ", ObarH1)

##### Subsection 2.2 Multiple runs of Ospats on sample
maxrun2 = 2*runmax
ObarS = zeros(Float64,1,maxrun2)
stratcy = Array{Int64,1}(n)
StratS =zeros(Int64,n,maxrun2)
cbObj = Array{Float64,1}(H)
cbObjS = zeros(Float64,H,maxrun2)

for run = 1:runmax
  TotTransf = 0
  println("       Start run nr  ", run)

##### Subsection 2.3. Initial stratification
  missing = n - H*floor(Int64, n/H)
  A = collect(1:H)
  B = vcat(A,A)
  repeat = n/H -2
  for rep = 1:repeat
    B = vcat(B,A)
  end
  fillup = collect(1:missing)
  B = vcat(B,fillup)
  strat0 = Array{Int64,1}(n)
  v = collect(1:n)
  w = v[randperm(rng,n)]
  for i = 1:n
    strat0[w[i]] = B[i]
  end

##### Subsection 2.4. Contributions from sample-strata
Sd2 = zeros(1,H)
for strat = 1:H
  Sd2[strat] = 0
  for i = 1:(n-1)
    if strat0[i] == strat
      for j = (i+1):n
        if strat0[j] == strat
          Sd2[strat] = Sd2[strat] + d2[i,j]
        end
      end
    end
  end
end
Sd2Init = Sd2
cbObj = sqrt(Sd2)
O = sum(cbObj)
ObarInit = O/n
# println("ObarInit = ", ObarInit)

##### Subsection 2.5. Transferring grid points
# println("Start transferring grid points")
stratcy = Array{Int64,2}(n,1)
stratcy = strat0
TotTransf = 0
TotCycle = 0
for cycle = 1:maxcycle
  transfers = 0
  u = randperm(n)
  for t = u
    Delta = 0
    change = 0
    A = stratcy[t]
    ij = find(stratcy .== A)
    dA = sum(d2[t,ij])
    sumd2tinA = dA
    Sd2Amint = Sd2[A] - sumd2tinA
    cbObjA = sqrt(abs(Sd2Amint))
    for stratnr = 1:H
      Delta = 0
      sumd2plus = 0
      if stratnr != A
        B = stratnr
        ij = find(stratcy .== B)
        dB = sum(d2[t,ij])
        sumd2plus = dB
        cbObjB = sqrt(abs(Sd2[B] + sumd2plus))
        Delta = cbObjA + cbObjB - cbObj[A] -cbObj[B]
        if Delta < O*1e-10
          change = 1
          transfers = transfers + 1
          stratcy[t] = B            # update stratification
          Sd2[A] = Sd2[A] - sumd2tinA
          Sd2[B] = Sd2[B] + sumd2plus
          cbObj = sqrt(abs(Sd2))
          Obj = sum(cbObj)
          Delta = 0
        end                       # if Delta < Obj*1e-10
      end                         # if stratnr != A
      if change ==1 break end
    end                           # for strat=1:H
  end                             # for t=u

  # println("cycle ", cycle, "     transfers = ", transfers)
  TotTransf = TotTransf + transfers
  if transfers == 0 break end     # stopping rule
  TotCycle = cycle
end                               # for cycle=1:maxcycle

# println("------------ End of transfer cycles ----------------")
println("Total number of transfers = ", TotTransf)
println("Number of iteration cycles = ", TotCycle)

O = sum(cbObj)
Obar_sample = O/n

##### Subsection 2.6. Save results from runs
ObarS = hcat(ObarS, Obar_sample)
StratS = hcat(StratS, stratcy)
cbObj = cbObj'
cbObjS = hcat(cbObjS, cbObj)

end                               # for run=1:runmax
# println("--------- End of runs----------")

first = 2*runmax +1
last = first + runmax -1
ObarSS = ObarS[first:last]
StratSS = StratS[1:n, first:last]
cbObjSS = cbObjS[1:H, first:last]
 println("Obar values from the runs: ",ObarSS)

##### Subsection 2.7 Select the best sample-stratification
best = sortperm(ObarSS)
strat_ordered = StratSS[1:n, best]
strat_best = strat_ordered[1:n, 1]
obar_ordered = ObarSS[best]
obar_best = obar_ordered[1]
cbObj_ordered = cbObjSS[1:H, best]
cbObj_best = cbObj_ordered[1:H, 1]

# println("cbObj_best = ", cbObj_best)

##### Subsection 2.8. Calculate size of sample-strata
Nh_sample = Array{Int64}(1,H)
for strat = 1:H
  k = find(strat_best .== strat)
  Nh_sample[strat] = length(k)
end

##### SECTION 3. ALLOCATION TO THE SAMPLE-STRATA
# println("-------------------- Start allocating ----")

##### Subsection 3.1 Define the initial grid stratification of
##### length N, with zeros for the non-sample points,
##### and strat_best for the sample points.
stratification = zeros(Int64,N,1)
stratification[sample] = strat_best

sumsample = cbObj_best.^2   # sums of distances within sample-strata

##### Subsection 3.2 Partition the sample data according to the
##### sample-stratification.
nh = Array{Int64}(1,H)
deltas = Array{Float64,2}(N,H)
deltas = zeros(N,H)

 for strat = 1:H
   k = find(stratification .== strat)
   nh[strat] = length(k)

##### Subsection 3.3 Extract a dataset for each variable from
##### sample-stratum "strat".
 xset = x[k]
 yset = y[k]
 z_predset = z_pred[k]
 s2set =s2[k]

##### Subsection 3.4 Calculate the sum of d2 for all gridpoints to
##### the points in sample-stratum "strat" using the datasets.
##### Calculate from the sum the increase of Obj (delta) if point
##### were allocated to "strat". Store deltas in N x H array "deltas".
for i =1:N
 sumd2 = 0
 delta = 0
 for j = 1:nh[strat]
    Lag = sqrt((x[i] - xset[j])^2 + (y[i] - yset[j])^2)
    Z2 = (z_pred[i] - z_predset[j])^2
    SS = s2[i] + s2set[j]
    Z2 = Z2/R2
    d2 = Z2 + SS - SS.*exp(-3*Lag/range)
    sumd2 = sumd2 + d2
 end    # j = 1:nh[strat]
delta = sqrt(sumsample[strat] + sumd2) - cbObj_best[strat]
deltas[i,strat] = delta

end   # for for i =1:N
end   # for strat = 1:H

##### Subsection 3.5 For each point t, find the stratum for which t
##### has the smallest delta, and allocate t to that stratum.
 for t = 1:N
 best = sortperm(deltas[t,1:H])
 stratification[t] = best[1]
 end    # for t = 1:N

 # println("---- End allocating ----")

 ##### SECTION 4. SAMPLE SIZES, CONTRIBUTIONS FROM GRID-STRATA TO O,
 ##### AND RELATIVE STANDARD ERROR

####                                 Size of the grid-strata:
Nh = Array{Int64}(1,H)
for strat = 1:H
  k = find(stratification .== strat)
  Nh[strat] = length(k)
end

### Contributions from grid-strata to O (as in Subsection 2.3)
# println("---- Calculating contribution from grid-strata to Obj ----")

Sd2 = zeros(1,H)
for i = 1:(N-1)
      for j = (i+1):N
        strati = stratification[i]
        stratj = stratification[j]
        if strati == stratj
          Lag = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
          Z2 = (z_pred[i] - z_pred[j])^2
          SS = s2[i] + s2[j]
          Z2 = Z2/R2
          d2 = Z2 + SS - SS.*exp(-3*Lag/range)
          Sd2[strati] = Sd2[strati] + d2
        end
      end
end
cbObj_grid = sqrt(Sd2)
O = sum(cbObj_grid)
Obar_grid = O/N

###                                            n_pred
meanzpred = sum(z_pred)/N
n_pred = (100*Obar_grid/(meanzpred*RSEmax))^2
n_pred = round(n_pred)
n_pred = convert(UInt64,n_pred)

###                                            Neyman allocation
sum_ahOh = 0
for h = 1:H
   ahOh = Nh[h]*cbObj_grid[h]
   sum_ahOh = sum_ahOh + ahOh
end

nh = Array{Float64}(1,H)
for h = 1:H
  nh[h] = n_pred *Nh[h] *cbObj_grid[h]/sum_ahOh
end
nh = round(nh)

###                                            Relative standard Error
RSE = 100*Obar_grid/(meanzpred*sqrt(nmax))

##### SECTION 5. STRATIFIED RANDOM SAMPLING WITH FINAL STRATIFICATION
##### AND NEYMAN ALLOCATION

n_tot = sum(nh)
n_tot = round(n_tot)
n_tot = convert(UInt64,n_tot)
points = Array{UInt64}(1,n_tot)
xs = Array{Float64}(1,n_tot)
ys = Array{Float64}(1,n_tot)

for h = 1:H
  k = find(stratification .== h)      # numbers of points in stratum h
  stratsize = length(k)
  v = collect(1:stratsize)
  w = v[randperm(rng, stratsize)]      # randomized indexes to points

  f=0
  for i = 1:nh[h]
    f=f+1                         # making Neyman allocations integer
  end

  k_rand = k[w]                   # randomized points in stratum h
  points_h = k_rand[1:f]          # put the first f points in sample

  if h == 1                       # concatenate all H vectors of points
    points = points_h
  elseif h > 1
    points = vcat(points,points_h)
  end
end
xs = x[points]                    # get x coordinates of sample points
ys = y[points]                    # get y coordinates of sample points
strata = stratification[points]       # get stratum numbers of sample points
sampnr = 1:n_tot                  # make sample numbers
strs = hcat(sampnr,strata,points,xs,ys)

##### SECTION 6. OUTPUT OF FINAL RESULTS

writedlm("Stratification", stratification)
writedlm("Sample", strs)

println("FINAL RESULTS: ")
println("Lowest Obar (O/n) from the sample = ", obar_best)
println("Obar (O/N) from the entire grid : ", Obar_grid)
println("Size of the sample-strata: ", Nh_sample)
println("Size of grid-strata : ", Nh)
println("Contributions to O from grid-strata (cbObj) : ", cbObj_grid)
println("Mean of z-predictions : ", meanzpred)
println("Predicted sample size needed to attain RSEmax : ", n_pred)
println("Neyman sample allocation : ", nh)
println("Predicted RSE for nmax : ", RSE)

println("---- END OF FUNCTION ALLOCATEMR ---- ")

end   # function allocatemr
