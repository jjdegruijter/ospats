
# This file "ospatsmr" contains the function ospatsmr().
# Version: 2017-10-28
# Author: Jaap de Gruijter

# CONTENT:
# 1. N x N matrix of generalized distances
# 2. Initial stratification
# 3. Contributions from strata to objective function (O)
# 4. Transfer grid points
# 5. Sample sizes and Relative Standard Error
# 6. Stratified random sampling
# 7. Output of final results

# For process monitoring, uncomment println() lines

function ospatsmr()
println("-------------------- START FUNCTION OSPATSMR ---")

##### SECTION 1. N x N MATRIX OF GENERALIZED DISTANCES

# println("----- Calculating matrix of generalized distances --")
Lags = zeros(N,N)
Z2 = zeros(N,N)
SS = zeros(N,N)
  for i = 1:(N-1)
  for j = (i+1):N
      Lags[i,j] = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
      Lags[j,i] = Lags[i,j]
      Z2[i,j] = (z_pred[i] - z_pred[j])^2
      Z2[j,i] = Z2[i,j]
      SS[i,j] = s2[i] + s2[j]
      SS[j,i] = SS[i,j]
  end
  end
  Z2 = Z2/R2
  # Cov = 0.5*SS.*exp(-3*Lags/ra)    explanation: 0.5*SS = Covmax
  # d2 = Z2 + SS - 2*Cov             explanation
d2 = Z2 + SS - SS.*exp(-3*Lags/range)

TOTd2 = sum(d2)/2
ObarH1 = sqrt(TOTd2)/N

##                                 Multiple runs of re-allocation
maxrun2 = 2*runmax
ObarS = zeros(Float64,1,maxrun2)
stratcy = Array{Int64,1}(N)
StratS = zeros(Int64,N,maxrun2)
cbObj = Array{Float64,1}(H)
cbObjS = zeros(Float64,H,maxrun2)

for run = 1:runmax
  TotTransf = 0
  println("       Start run nr  ", run)

##### SECTION 2. INITIAL STRATIFICATION

  # println("----- Generating random initial stratification ---")
  missing = N - H*floor(Int64, N/H)
  A = collect(1:H)
  B = vcat(A,A)
  repeat = N/H -2
  for rep = 1:repeat
    B = vcat(B,A)
  end
  fillup = collect(1:missing)
  B = vcat(B,fillup)
  strat0 = Array{Int64,1}(N)
  v = collect(1:N)
  w = v[randperm(rng, N)]
  for i = 1:N
    strat0[w[i]] = B[i]
  end

##### SECTION 3. CONTRIBUTIONS FROM STRATA TO OBJECTIVE FUNCTION

# println("--- Calculating contributions from strata to O ---")
Sd2 = zeros(1,H)
for strat = 1:H
  Sd2[strat] = 0
  for i = 1:(N-1)
    if strat0[i] == strat
      for j = (i+1):N
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
ObarInit = O/N
# println("ObarInit = ", ObarInit)

##### SECTION 4. TRANSFER GRID POINTS

# println("--------------- Transferring grid points ---")
stratcy = Array{Int64,2}(N,1)
stratcy = strat0
TotTransf = 0
TotCycle = 0
for cycle = 1:maxcycle
  transfers = 0
  u = randperm(N)
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
          O = sum(cbObj)
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

println("Total number of transfers = ", TotTransf)
println("Number of iteration cycles = ", TotCycle)

O = sum(cbObj)
ObarFinal = O/N

##                                       save results from runs
ObarS = hcat(ObarS, ObarFinal)
StratS = hcat(StratS, stratcy)
cbObj = cbObj'
cbObjS = hcat(cbObjS, cbObj)

end                               # for run=1:runmax
# println("--------- end of runs----------")

first = 2*runmax +1
last = first + runmax -1
ObarSS = ObarS[first:last]
StratSS = StratS[1:N, first:last]
cbObjSS = cbObjS[1:H, first:last]

##                                      select the best stratification
best = sortperm(ObarSS)
strat_ordered = StratSS[1:N, best]
strat_best = strat_ordered[1:N, 1]
cbObj_ordered = cbObjSS[1:H, best]
cbObj_best = cbObj_ordered[1:H, 1]
obar_ordered = ObarSS[best]
obar_best = obar_ordered[1]

##### SECTION 5. SAMPLE SIZES AND RELATIVE STANDARD ERROR

Nh = Array{Int64}(1,H)
for h = 1:H
  k = find(strat_best .== h)
  Nh[h] = length(k)
end

###                                            n_pred
meanzpred = sum(z_pred)/N
n_pred = (100*obar_best/(meanzpred*RSEmax))^2
n_pred = round(n_pred)
n_pred = convert(UInt64,n_pred)

sum_ahOh = 0
for h = 1:H
   ahOh = Nh[h]*cbObj_best[h]
   sum_ahOh = sum_ahOh + ahOh
end

###                                            Neyman allocation
nh = Array{Float64}(1,H)
for h = 1:H
  nh[h] = n_pred *Nh[h] *cbObj_best[h]/sum_ahOh
end
nh = round(nh)

###                                            Relative standard Error
RSE = 100*obar_best/(meanzpred*sqrt(nmax))

##### SECTION 6. STRATIFIED RANDOM SAMPLING WITH FINAL STRATIFICATION
##### AND NEYMAN ALLOCATION

n_tot = sum(nh)
n_tot = round(n_tot)
n_tot = convert(UInt64,n_tot)
points = Array{UInt64}(1,n_tot)
xs = Array{Float64}(1,n_tot)
ys = Array{Float64}(1,n_tot)

for h = 1:H
  k = find(strat_best .== h)      # numbers of points in stratum h
  stratsize = length(k)
  v = collect(1:stratsize)
  w = v[randperm(rng,stratsize)]  # randomized indexes to points

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
strata = strat_best[points]       # get stratum numbers of sample points
sampnr = 1:n_tot                  # make sample numbers
strs = hcat(sampnr,strata,points,xs,ys)

##### SECTION 7. OUTPUT OF FINAL RESULTS

writedlm("Stratification", strat_best)
writedlm("Sample", strs)

println("FINAL RESULTS: ")
println("Value of Obar (O/N) without stratification : ", ObarH1)
println("Values of Obar : ", ObarSS)
println("Lowest value of Obar : ", obar_best)
println("Sizes of grid-strata : ", Nh)
println("Contributions to O from strata (cbObj) : ", cbObj_best)
println("Mean of z-predictions : ", meanzpred)
println("Predicted sample size needed to attain RSEmax : ", n_pred)
println("Neyman sample allocation : ", nh)
println("Predicted RSE for nmax : ", RSE)

println("--------------------- END FUNCTION OSPATSMR ---")
end                         # function ospatsmr
