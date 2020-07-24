## Importing the Required Libraries
import matplotlib.pyplot as plt
import math
from math import pow
import numpy
import numpy as np
from statistics import harmonic_mean, geometric_mean
import operator
import csv
import pandas as pd

## Creating the global functions required to read and convert the data
#Reverse Generator from MBVD Parameters
def genMBVDfrq(Rs, Rm, Lm, Cm, R0, C0, y_abs_gen, y_real_gen, y_imag_gen, z_abs_gen, z_real_gen, z_imag_gen, frq_gen):
  ## Static Equations for Q and kt
  kt = ((pow((math.pi), 2) / 8) * (Cm / C0))
  Q = ((1 / Rm) * (pow((Lm / Cm), 0.5)))

  ## Outputing Q and kt
  print("Generated Results:")
  print(" ")
  print("Kt-sqr is " + str(kt))
  print("Q is " + str(Q))
  print(" ")

  freq = 150000000
  #x_start = 150000000
  freq_end = 165000000
  x = []

  ## For each input, we generate the frequencies at the correct hz and perform calculations to find other variables
  while freq <= freq_end:
    x.append(freq)
    ## Math to find the other variables based on pdf
    p = (R0 + Rm)
    a0 = (1.0 / (2 * math.pi * freq * C0))
    am = (2 * math.pi * freq * Lm) - (1.0 / (2 * math.pi * freq * Cm))
    W = (am - a0)
    Xt = ((Rs * R0) + (Rs * Rm) + (R0 * Rm) + (a0 * am))
    Yt = ((Rs * am) - (Rs * a0) + (R0 * am) - (Rm * a0))
    mt = (1.0 / ((Xt * Xt) + (Yt * Yt)))
    admitance_r = ((p * Xt) + (W * Yt)) * mt
    admitance_i = ((W * Xt) - (p * Yt)) * mt
    zt = (1.0 / ((p * p) + (W * W)))
    impedance_r = zt * admitance_r * (1.0 / mt)
    impedance_i = -1.0 * zt * admitance_i * (1.0 / mt)

    admitance_abs = math.sqrt((admitance_r * admitance_r) + (admitance_i * admitance_i))

    ##Appending the values to the list in order to keep running tabs
    y_abs_gen.append(admitance_abs)
    y_real_gen.append(admitance_r)
    y_imag_gen.append(admitance_i)
    z_abs_x = math.sqrt((impedance_r * impedance_r) + (impedance_i * impedance_i))
    z_abs_gen.append(z_abs_x)
    z_real_gen.append(impedance_r)
    z_imag_gen.append(impedance_i)
    frq_gen.append(freq)
    freq = freq + 937.5
  return x

#Reverse Generator from MBVD Parameters
def genMBVD(Rs, Rm, Lm, Cm, R0, C0, y_abs_gen, y_real_gen, y_imag_gen, z_abs_gen, z_real_gen, z_imag_gen, frq_gen):
  ## Static Equations for Q and kt
  kt = ((pow((math.pi), 2) / 8) * (Cm / C0))
  Q = ((1 / Rm) * (pow((Lm / Cm), 0.5)))

  ## Outputing Q and kt
  print("Generated Results:")
  print(" ")
  print("Kt-sqr is " + str(kt))
  print("Q is " + str(Q))
  print(" ")

  x = [n for n in range(150000, 164000)]

  ## For each input, we generate the frequencies at the correct hz and perform calculations to find other variables
  for num in x:
    # Input Frequencies
    freq = num * pow(10, 3)

    ## Math to find the other variables based on pdf
    p = (R0 + Rm)
    a0 = (1.0 / (2 * math.pi * freq * C0))
    am = (2 * math.pi * freq * Lm) - (1.0 / (2 * math.pi * freq * Cm))
    W = (am - a0)
    Xt = ((Rs * R0) + (Rs * Rm) + (R0 * Rm) + (a0 * am))
    Yt = ((Rs * am) - (Rs * a0) + (R0 * am) - (Rm * a0))
    mt = (1.0 / ((Xt * Xt) + (Yt * Yt)))
    admitance_r = ((p * Xt) + (W * Yt)) * mt
    admitance_i = ((W * Xt) - (p * Yt)) * mt
    zt = (1.0 / ((p * p) + (W * W)))
    impedance_r = zt * admitance_r * (1.0 / mt)
    impedance_i = -1.0 * zt * admitance_i * (1.0 / mt)

    admitance_abs = math.sqrt((admitance_r * admitance_r) + (admitance_i * admitance_i))

    ##Appending the values to the list in order to keep running tabs
    y_abs_gen.append(admitance_abs)
    y_real_gen.append(admitance_r)
    y_imag_gen.append(admitance_i)
    z_abs_x = math.sqrt((impedance_r * impedance_r) + (impedance_i * impedance_i))
    z_abs_gen.append(z_abs_x)
    z_real_gen.append(impedance_r)
    z_imag_gen.append(impedance_i)
    frq_gen.append(freq)
  return x


##Plotting the Data
def matplot(title, x, y, x_label, y_label, number, multi_x=None, multi_y=None, label1=None, label2=None):
  plt.title(title)
  if not (multi_x == None or multi_y == None):
    plt.plot(x, y, label = str(label1))
    plt.plot(multi_x, multi_y, label = str(label2))
    plt.legend()
  else:
    plt.plot(x, y)
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  plt.grid(True)
  #plt.show()
  plt.savefig("" + str(title) + "_" + str(number) + ".png")
  plt.close()
  
## Finding Resistor Values
def resTest(max_y_ind, max_z_ind, max_y_val, max_z_val, y_mean, z_mean, y_abs, z_abs, ws, wp, Cm_estimate, sum_of_rs_and_r0_est, ydb1, ydb2, zdb1, zdb2):
  y_threshold = ((max_y_val - y_mean) * math.sqrt(.5)) + y_mean
  z_threshold = ((max_z_val - z_mean) * math.sqrt(.5)) + z_mean

  print("Y Thres:" + str(y_threshold))
  print("Z Thres:" + str(z_threshold))

  p = o = max_y_ind
  for num in range(max_y_ind, len(y_abs)):
      temp_ind = num - max_y_ind
      if (y_abs[num] >= y_threshold):
          p = num
      if (y_abs[max_y_ind - temp_ind] >= y_threshold):
          o = max_y_ind - temp_ind

  y_width1 = frq[p] - frq[o]

  q = w = max_z_ind
  e = 0
  for num in range(max_z_ind, len(z_abs)):
    temp_in = num - max_z_ind
    if (z_abs[num] >= z_threshold):
      q = num
    if (z_abs[max_z_ind - temp_in] >= z_threshold):
        w = max_z_ind - temp_in
        if (abs(e) == temp_in):
          e = e - 1
          w = max_z_ind - temp_in

  z_width1 = frq[q] - frq[w]

  print("\nY Frq1:" + str(frq[p]))
  print("Y Frq2:" + str(frq[o]))

  print("\nY Frq1 Diff:" + str((frq[p] - frq[ydb2])))
  print("Y Frq2 Diff:" + str((frq[o] - frq[ydb1])))

  print("\nY Val1:" + str(y_abs[p]))
  print("Y Val2:" + str(y_abs[o]))

  print("\nY Val1 Diff:" + str((y_abs[p] - y_abs[ydb2])))
  print("Y Val2 Diff:" + str((y_abs[o] - y_abs[ydb1])))

  print("\nZ Frq1:" + str(frq[q]))
  print("Z Frq2:" + str(frq[w]))

  print("\nZ Frq1 Diff:" + str((frq[q] - frq[zdb2])))
  print("Z Frq2 Diff:" + str((frq[w] - frq[zdb1])))

  print("\nZ Val1:" + str(z_abs[q]))
  print("Z Val2:" + str(z_abs[w]))

  print("\nZ Val1 Diff:" + str((z_abs[q] - z_abs[zdb2])))
  print("Z Val2 Diff:" + str((z_abs[w] - z_abs[zdb1])))

  print("\nY Width:" + str(y_width1))
  print("Z Width:" + str(z_width1))

  Qs0_estimate1 = ws / (y_width1 * 2 * math.pi)
  Qp0_estimate1 = wp / (z_width1 * 2 * math.pi)

  print("Qs0: " + str(Qs0_estimate1))
  print("Qp0: " + str(Qp0_estimate1))

  rm_and_rs_est = 1.0 / (ws * Cm_estimate * Qs0_estimate1)
  rm_and_r0_est = 1.0 / (wp * Cm_estimate * Qp0_estimate1)

  A = numpy.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
  B = numpy.array([rm_and_rs_est, sum_of_rs_and_r0_est, rm_and_r0_est])
  R_list = numpy.linalg.solve(A,B)

  r0_est1 = R_list[2]
  rm_est1 = R_list[0]
  rs_est1 = R_list[1]
  
  return(rm_and_rs_est, rm_and_r0_est, r0_est1, rm_est1, rs_est1, Qs0_estimate1, Qp0_estimate1)

## Converting cartesian to polar
def cart2pol(x, y):
  rho = numpy.sqrt(x**2 + y**2)
  phi = numpy.arctan2(y, x)
  return numpy.array((rho, phi))

## Converting polar to cartesian
def pol2cart(rho, phi):
  x = rho * numpy.cos(phi)
  y = rho * numpy.sin(phi)
  return numpy.array((x, y))

#Converting Cartesian inputs into magnitude
def c2m(x, y, d, c):
  top = np.sqrt((x*x) + (y*y))
  num = c * top
  total = num/d
  return total

#Converting db back to standard numbers
def db2lin(A, P1 = 1):
  x = math.pow(10,(A/20)) # Used to be A/10
  P2 = P1 * x
  return P2

#Converting amplitude to db
def lin2db(A, P1 = 1):
  db = 20 * math.log10((A/P1))
  return db

#Converting amplitude to db for lists
def lin2db_list(inlist):
  outlist = []
  for num in inlist:
    outlist.append(lin2db(num))
  return outlist

#Converting Polar Inputs into Cartesian coordinates
def p2c(rho, phi):
    phi = np.radians(phi)
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y

data = pd.read_csv(".//input.csv")
x = data['Freq(Hz)'].tolist()

#Creating a Dataframe to Save File
dataframe = pd.DataFrame()

# Loading the data into arrays
S11db = data['S11(DB)'].tolist()
S11deg = data['S11(DEG)'].tolist()
S12db = data['S12(DB)'].tolist()
S12deg = data['S12(DEG)'].tolist()
S21db = data['S21(DB)'].tolist()
S21deg = data['S21(DEG)'].tolist()
S22db = data['S22(DB)'].tolist()
S22deg = data['S22(DEG)'].tolist()
EQdb = data['Eq=(1-S11)/(1+S11)(DB)'].tolist()
EQdeg = data['Eq=(1-S11)/(1+S11)(DEG)'].tolist()

y_eq_abs = []
y_eq_real = []
y_eq_imag = []
y_abs = []
y_real = []
y_imag = []
z_abs = []
z_real = []
z_imag = []
frq = []
positive = []
zerox = []

#Math based on the information in the Google Docs
for index, item in enumerate(S11db):  
  S11_mag = db2lin(S11db[index])
  
  A, B = p2c(S11_mag, S11deg[index])

  eq_mag = db2lin(EQdb[index])

  C, D = p2c(eq_mag, EQdeg[index])

  paper_eq1 = (1-A*A-B*B)
  paper_eq2 = (-2*B)
  paper_eqden = (A*A + 2*A + 1 + B*B)

  y_eq = math.sqrt((C*C)+(D*D))

  z_c = 1 - (2*A) + (A*A) + (B*B)
  z_r = (50 * (1 - (A*A) - (B*B)))/z_c
  z_i = (50 * (2*B))/z_c
  
  z_mag = math.sqrt((z_r*z_r) + (z_i*z_i))

  Z_test = (1/50)

  Y_mag = c2m(paper_eq1, paper_eq2, paper_eqden, Z_test)
 
  #Adding each admitance to list y
  y_eq_abs.append(y_eq)
  y_eq_real.append(C)
  y_eq_imag.append(D)
  y_abs.append(Y_mag)
  z_abs.append(z_mag)
  y_real.append(Z_test*paper_eq1/paper_eqden)
  y_imag.append(Z_test*paper_eq2/paper_eqden)
  frq.append(x[index])
  z_real.append(z_r)
  z_imag.append(z_i)
  positive.append((z_i)>0)

## The designated off frequency values that are used as reference
off_freq = [1, 4800, 16000]

print("Off-Freq: " + str(frq[off_freq[0]]) + ", " + str(frq[off_freq[1]]) + ", "+ str(frq[off_freq[2]]))
print("Off-Freq Z Real Value: " + str(z_real[off_freq[0]]) + ", " + str(z_real[off_freq[1]]) + ", "+ str(z_real[off_freq[2]]))
## In the above list, these are the indcies of the off frequency values

##Setting Up to Write to CSV
labels = ['Frequency', 'admitance_abs', 'admitance_r', 'admitance_i', 'impedance_abs', 'impedance_r', 'impedance_i']

rows = zip(frq, y_abs, y_real, y_imag, z_abs, z_real, z_imag)

newfilePath = "data_converted.csv"

##Writing to a File
with open(newfilePath, "w") as f:
    writer = csv.writer(f)
    writer.writerow(labels)
    for row in rows:
        writer.writerow(row)

##Plotting the Data
plotcount = 0
matplot('Impedance vs Frequency', frq, z_abs, 'Frequency', 'Impedance', plotcount)
matplot('Admittance vs Frequency', frq, y_abs, 'Frequency', 'Admittance', plotcount)

matplot('Admittance vs Frequency (EQ)', frq, y_eq_abs, 'Frequency', 'Admittance', plotcount)

matplot('Impedance (Imaginary) vs Frequency', frq, z_imag, 'Frequency', 'Impedance', plotcount)
matplot('Impedance (Real) vs Frequency', frq, z_real, 'Frequency', 'Impedance', plotcount)

matplot('Admittance (Imaginary) vs Frequency', frq, y_imag, 'Frequency', 'Admittance', plotcount)
matplot('Admittance (Real) vs Frequency', frq, y_real, 'Frequency', 'Admittance', plotcount)


y_abs_db = lin2db_list(y_abs)
z_abs_db = lin2db_list(z_abs)

matplot('Impedance (Gen DB) vs Frequency', frq, z_abs_db, 'Frequency', 'Impedance', plotcount)
matplot('Admittance (Gen DB) vs Frequency', frq, y_abs_db, 'Frequency', 'Admittance', plotcount)

plotcount = plotcount + 1

max_ydb_ind, max_ydb_val = max(enumerate(y_abs_db), key=operator.itemgetter(1))
max_zdb_ind, max_zdb_val = max(enumerate(z_abs_db), key=operator.itemgetter(1))

y_abs_db = [-i for i in y_abs_db]

ydb_mean = -1 * geometric_mean(y_abs_db)
zdb_mean = geometric_mean(z_abs_db)

y_abs_db = [-i for i in y_abs_db]

frq_y_db_1 = y_abs_db.index(min(y_abs_db[:max_ydb_ind], key=lambda x:abs(x-(max_ydb_val - 3))))
frq_y_db_2 = y_abs_db.index(min(y_abs_db[max_ydb_ind:], key=lambda x:abs(x-(max_ydb_val - 3))))
frq_z_db_1 = z_abs_db.index(min(z_abs_db[:max_zdb_ind], key=lambda x:abs(x-(max_zdb_val - 3))))
frq_z_db_2 = z_abs_db.index(min(z_abs_db[max_zdb_ind:], key=lambda x:abs(x-(max_zdb_val - 3))))

print("Y DB FRQ #1 Comp: " + str(frq[frq_y_db_1]))
print("Y DB FRQ #2 Comp: " + str(frq[frq_y_db_2]))
print("Z DB FRQ #1 Comp: " + str(frq[frq_z_db_1]))
print("Z DB FRQ #2 Comp: " + str(frq[frq_z_db_2]))

## Let's estimate the values from the data created
## Reverse engieering the values from the generated data

## Finding the C0 estimate
## Formula: l/avg(off_freq)
sum_of_imag_at_off_freq = -2 * math.pi * (
    (frq[int(off_freq[0])] * (z_imag[int(off_freq[0])])) +
    (frq[int(off_freq[1])] * (z_imag[int(off_freq[1])])) +
    (frq[int(off_freq[2])] * (z_imag[int(off_freq[2])])))

sum_of_imag_at_off_freq = sum_of_imag_at_off_freq / 3.0
C0_estimate = 1.0 / (sum_of_imag_at_off_freq)

## Outputing the results
print("C0_estimate: " + str(C0_estimate))
print(" ")

## Finding the sum of R0 and Rs
sum_of_rs_and_r0_est = (abs(z_real[int(off_freq[0])]) + abs(z_real[int(
    off_freq[1])]) + abs(z_real[int(off_freq[2])]))

sum_of_rs_and_r0_est = sum_of_rs_and_r0_est / 3.0

zerox = numpy.where(numpy.bitwise_xor(positive[1:], positive[:-1]))[0]
fs = frq[zerox[0]]
fp = frq[zerox[1]]
print("Fs: " + str(fs))
print("Fp: " + str(fp))
ws = 2 * fs * math.pi
wp = 2 * fp * math.pi
r = 1.0 / (((fp / fs)**2) - 1)

Cm_estimate = C0_estimate / r
print("Cm_estimate: " + str(Cm_estimate))
print(" ")
if len(zerox) > 2:
  print(
    "More than 2 zero crossings reported, need to re evaluate the zero crossings"
    )
  print(" ")

Lm_estimate = 1.0 / (ws**2 * Cm_estimate)
print("Lm_estimate: " + str(Lm_estimate))
print(" ")

#Finding Fs and bandwitdh measurements
max_y_ind, max_y_val = max(enumerate(y_abs), key=operator.itemgetter(1))
max_z_ind, max_z_val = max(enumerate(z_abs), key=operator.itemgetter(1))
min_y_ind, min_y_val = min(enumerate(y_abs), key=operator.itemgetter(1))
min_z_ind, min_z_val = min(enumerate(z_abs), key=operator.itemgetter(1))

#y_mean = harmonic_mean(y_abs)
#z_mean = harmonic_mean(z_abs)

y_mean = geometric_mean(y_abs)
z_mean = geometric_mean(z_abs)

rm_and_rs_est, rm_and_r0_est, r0_est, rm_est, rs_est, QS0, QP0 = resTest(max_y_ind, max_z_ind, max_y_val, max_z_val, y_mean, z_mean, y_abs, z_abs, ws, wp, Cm_estimate, sum_of_rs_and_r0_est, frq_y_db_1, frq_y_db_2, frq_z_db_1, frq_z_db_2)

print("__________________________________________________________________ ")
print("Further Testing - 50% reduction from peak bandwidth test")
print("sum_rm_and_rs_est: " + str(rm_and_rs_est))
print("sum_rm_and_r0_est: " + str(rm_and_r0_est))
print("sum_rs_and_r0_est: " + str(sum_of_rs_and_r0_est))
print(" ")
print("r0_est: " + str(r0_est))
print("rm_est: " + str(rm_est))
print("rs_est: " + str(rs_est))
print(" ")
print("Mean Approximation Y: " + str(y_mean))
print("Mean Approximation Z: " + str(z_mean))
print("Qs0: " + str(QS0))
print("Qp0: " + str(QP0))
print(" ")
print("__________________________________________________________________ ")

print("More refinement is still possible...\n\n ")

## Instantiating lists to hold generated values
y_abs_gen = []
y_real_gen = []
y_imag_gen = []
z_abs_gen = []
z_real_gen = []
z_imag_gen = []
frq_gen = []

Rs = rs_est
Rm = rm_est
Lm = Lm_estimate
Cm = Cm_estimate
R0 = r0_est
C0 = C0_estimate

x = genMBVD(Rs, Rm, Lm, Cm, R0, C0, y_abs_gen, y_real_gen, y_imag_gen, z_abs_gen, z_real_gen, z_imag_gen, frq_gen)

##Plotting the Data
matplot('Impedance vs Frequency (Generated)', x, z_abs_gen, 'Frequency (KHz)', 'Impedance', plotcount)
matplot('Admittance vs Frequency (Generated)', x, y_abs_gen, 'Frequency (KHz)', 'Admittance', plotcount)

x = [i * 1000 for i in x]

matplot('Impedance vs Frequency [Abs] (Multi)', x, z_abs_gen, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_abs, 'Generated', 'Original')
matplot('Admittance vs Frequency [Abs] (Multi)', x, y_abs_gen, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_abs, 'Generated', 'Original')

matplot('Impedance vs Frequency [Real] (Multi)', x, z_real_gen, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_real, 'Generated', 'Original')
matplot('Admittance vs Frequency [Real] (Multi)', x, y_real_gen, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_real, 'Generated', 'Original')

matplot('Impedance vs Frequency [Img] (Multi)', x, z_imag_gen, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_imag, 'Generated', 'Original')
matplot('Admittance vs Frequency [Img] (Multi)', x, y_imag_gen, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_imag, 'Generated', 'Original')

##Setting Up to Write to CSV
labels = ['Frequency', 'admitance_abs', 'admitance_r', 'admitance_i', 'impedance_abs', 'impedance_r', 'impedance_i']

rows = zip(x, y_abs_gen, y_real_gen, y_imag_gen, z_abs_gen, z_real_gen, z_imag_gen)

newfilePath = "data_gen_converted.csv"

##Writing to a File
with open(newfilePath, "w") as f:
    writer = csv.writer(f)
    writer.writerow(labels)
    for row in rows:
        writer.writerow(row)

y_abs_genfrq = []
y_real_genfrq = []
y_imag_genfrq = []
z_abs_genfrq = []
z_real_genfrq = []
z_imag_genfrq = []
frq_genfrq = []

###Something here: we need to edit below
xfrq = genMBVDfrq(Rs, Rm, Lm, Cm, R0, C0, y_abs_genfrq, y_real_genfrq, y_imag_genfrq, z_abs_genfrq, z_real_genfrq, z_imag_genfrq, frq_genfrq)

matplot('Impedance vs Frequency [Abs] (Multi-FRQadj)', xfrq, z_abs_genfrq, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_abs, 'Generated', 'Original')
matplot('Admittance vs Frequency [Abs] (Multi-FRQadj)', xfrq, y_abs_genfrq, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_abs, 'Generated', 'Original')

matplot('Impedance vs Frequency [Real] (Multi-FRQadj)', xfrq, z_real_genfrq, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_real, 'Generated', 'Original')
matplot('Admittance vs Frequency [Real] (Multi-FRQadj)', xfrq, y_real_genfrq, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_real, 'Generated', 'Original')

matplot('Impedance vs Frequency [Img] (Multi-FRQadj)', xfrq, z_imag_genfrq, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_imag, 'Generated', 'Original')
matplot('Admittance vs Frequency [Img] (Multi-FRQadj)', xfrq, y_imag_genfrq, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_imag, 'Generated', 'Original')

##Setting Up to Write to CSV
labels = ['Frequency', 'admitance_abs', 'admitance_r', 'admitance_i', 'impedance_abs', 'impedance_r', 'impedance_i']

rows = zip(xfrq, y_abs_genfrq, y_real_genfrq, y_imag_genfrq, z_abs_genfrq, z_real_genfrq, z_imag_genfrq)

newfilePath = "data_gen_converted_frq_adj.csv"

##Writing to a File
with open(newfilePath, "w") as f:
    writer = csv.writer(f)
    writer.writerow(labels)
    for row in rows:
        writer.writerow(row)
###End the something here


plotcount = plotcount + 1

s = 0
p = 0

cur_qs = QS0 + s
cur_qp = QP0 + p 

print("Current QS0: " + str(cur_qs))
print("Current QP0: " + str(cur_qp))

rate = 0.01
#.00000075
#.0000003
#.00000025 - Okay - time taken
#.000001 - too fast
#.0000001 # Learning rate - working but slow

precision = 3000
iters = 0
cost_f = 0
cost_ds = 0
cost_dp = 0
temp = 0
xp = 1/(wp * C0)

while ((cost_f > precision) or (temp == 0)):
  temp += 1
  for ind, frequ in enumerate(frq):
    S = frequ/ws
    P = frequ/wp
    iters += 1
    
    ds = ((-1 * xp * cur_qp) / (P * (cur_qs**2))) * ((cur_qp * (1 - (P**2)) * S) / (((cur_qp**2) * ((1 - (P**2))**2)) + (P**2)))
    dp = (((xp)/(P * cur_qs))*(((2 * cur_qp * (1 - (P**2)) * S) - (cur_qs * P * (1 - (S**2))))/(((cur_qp**2)*((1-(P**2))**2)) + (P**2)))) + (((2 * xp * (cur_qp**2) * ((1 - (P**2))**2))/(P * cur_qs))*(((cur_qp * (1 - (P**2)) * S) - (cur_qs * P * (1 - (S**2))))/((((cur_qp**2)*((1-(P**2))**2)) + (P**2))**2)))
    
    ind_freq = frq_gen.index(min(frq_gen, key=lambda x:abs(x-frequ)))

    cost_f += (z_real[ind] - z_real_gen[ind_freq])**2
    cost_ds += (z_real[ind] - z_real_gen[ind_freq])*(ds)
    cost_dp += (z_real[ind] - z_real_gen[ind_freq])*(dp)

  cost_f = cost_f/(2*(iters))
  if temp == 1:
    print("Intial Cost: "+ str(cost_f)) 
  gradient_cost_ds = -1 * cost_ds/(iters)
  gradient_cost_dp = -1 * cost_dp/(iters)

  delta_cost_ds = -rate * (gradient_cost_ds)**2
  delta_cost_dp = -rate * (gradient_cost_dp)**2

  delta_cost = delta_cost_ds + delta_cost_dp

  cost_f = cost_f + delta_cost
  cur_qs = cur_qs - (rate * (gradient_cost_ds))
  cur_qp = cur_qp - (rate * (gradient_cost_dp))
  
  print("\nCost: " + str(cost_f))
  print("QS0: " + str(cur_qs))
  print("QP0: " + str(cur_qp))

  iters = 0
  rm_and_rs_est = 1.0 / (ws * Cm * cur_qs)
  rm_and_r0_est = 1.0 / (wp * Cm * cur_qp)

  A = numpy.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
  B = numpy.array([rm_and_rs_est, sum_of_rs_and_r0_est, rm_and_r0_est])
  R_list = numpy.linalg.solve(A,B)

  R0 = R_list[2]
  Rm = R_list[0]
  Rs = R_list[1]

  if temp == 1:
    print("First Rs:" + str(Rs))
    print("First Rm:" + str(Rm))
    print("First R0:" + str(R0))

  x = [n for n in range(150000, 164000)]

  ## Instantiating lists to hold generated values
  z_real_gen = []
  frq_gen = []

  ## For each input, we generate the frequencies at the correct hz and perform calculations to find other variables
  for num in x:
    # Input Frequencies
    freq = num * pow(10, 3)

    ## Math to find the other variables based on pdf
    p = (R0 + Rm)
    a0 = (1.0 / (2 * math.pi * freq * C0))
    am = (2 * math.pi * freq * Lm) - (1.0 / (2 * math.pi * freq * Cm))
    W = (am - a0)
    Xt = ((Rs * R0) + (Rs * Rm) + (R0 * Rm) + (a0 * am))
    Yt = ((Rs * am) - (Rs * a0) + (R0 * am) - (Rm * a0))
    mt = (1.0 / ((Xt * Xt) + (Yt * Yt)))
    admitance_r = ((p * Xt) + (W * Yt)) * mt
    zt = (1.0 / ((p * p) + (W * W)))
    impedance_r = zt * admitance_r * (1.0 / mt)
    z_real_gen.append(impedance_r)
    frq_gen.append(freq)

print("\nFinal Cost: " + str(cost_f))
print("Final QS0: " + str(cur_qs))
print("Final QP0: " + str(cur_qp))
print("Final Rs:" + str(Rs))
print("Final Rm:" + str(Rm))
print("Final R0:" + str(R0))

## Instantiating lists to hold generated values
y_abs_gen2 = []
z_abs_gen2 = []
frq_gen2 = []
y_real_gen2 = []
y_imag_gen2 = []
z_real_gen2 = []
z_imag_gen2 = []

x = genMBVD(Rs, Rm, Lm, Cm, R0, C0, y_abs_gen2, y_real_gen2, y_imag_gen2, z_abs_gen2, z_real_gen2, z_imag_gen2, frq_gen2)

x = [i * 1000 for i in x]

matplot('Impedance vs Frequency (Multi)', x, z_abs_gen2, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_abs, 'Generated', 'Original')
matplot('Admittance vs Frequency (Multi)', x, y_abs_gen2, 'Frequency (KHz)', 'Admittance', plotcount, frq, y_abs, 'Generated', 'Original')
matplot('Impedance (Real) vs Frequency (Multi)', x, z_real_gen, 'Frequency (KHz)', 'Impedance', plotcount, frq, z_real, 'Generated', 'Original')

print("Done")
