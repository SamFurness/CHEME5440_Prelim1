#Sam Furness
#CHEME5440 Prelim1 Question 3 Part A,B

include("Flux.jl")
#We are passing arrays of all the constraints lower and upper bounds
  #as well as the stoichiometric array in this FBA function
using LinearAlgebra
using Conda
using PyPlot

#Define the stoichiometrix matrix. 17species x 15fluxes. Note that the transcript is 924nt
  #and that the peptide is 308aa. So n=924 and a=308 in the Allen/Pallson Reactions
Stoichiometric_matrix = Array{Float64,2}(undef,17,15)
Stoichiometric_matrix[1,:]=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0];#gene G
Stoichiometric_matrix[2,:]=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0];#RNAP
Stoichiometric_matrix[3,:]=[1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0];#gene+RNAP G*
Stoichiometric_matrix[4,:]=[0 -924 0 0 0 0 0 1 0 0 0 0 0 0 0];#NTP
Stoichiometric_matrix[5,:]=[0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0];#mRNA
Stoichiometric_matrix[6,:]=[0 2*924 0 0 2*308 2 0 0 0 0 0 0 0 0 -1];#Pi
Stoichiometric_matrix[7,:]=[0 0 924 0 0 0 0 0 0 -1 0 0 0 0 0];#NMP
Stoichiometric_matrix[8,:]=[0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0];#ribosome rib
Stoichiometric_matrix[9,:]=[0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0];#ribosome+mRNA rib*
Stoichiometric_matrix[10,:]=[0 0 0 0 -308 1 0 0 0 0 0 0 0 0 0];#AAtRNA
Stoichiometric_matrix[11,:]=[0 0 0 0 -2*308 0 0 0 0 0 0 0 1 0 0];#GTP
Stoichiometric_matrix[12,:]=[0 0 0 0 308 -1 0 0 0 0 0 0 0 0 0];#tRNA
Stoichiometric_matrix[13,:]=[0 0 0 0 2*308 0 0 0 0 0 0 0 0 -1 0];#GDP
Stoichiometric_matrix[14,:]=[0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0];#protein
Stoichiometric_matrix[15,:]=[0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0];#AA
Stoichiometric_matrix[16,:]=[0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0];#ATP
Stoichiometric_matrix[17,:]=[0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0];#AMP
#Stoichiometric_matrix[18,:]=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];#AA_external
#Stoichiometric_matrix[19,:]=[0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];#NTP_external
#Stoichiometric_matrix[20,:]=[0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0];#protein_external
#Stoichiometric_matrix[21,:]=[0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0];#NMP_external
#Stoichiometric_matrix[22,:]=[0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];#ATP_external
#Stoichiometric_matrix[23,:]=[0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0];#AMP_external
#Stoichiometric_matrix[24,:]=[0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0];#GTP_external
#Stoichiometric_matrix[25,:]=[0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0];#GDP_external
#Stoichiometric_matrix[26,:]=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1];#Pi_external

#Writing an expression for transcription and translation
  #Define parameters
  eX=60.0;#transcription elongation rate nt/s
  eL=16.5#translation elongation rate aa/s
  LX=924.0;#length of mRNA transcript nt
  LL=308.0;#length of peptide aa
  kEX=eX/LX;#s^-1
  kEL=eL/LL;#s^-1
  Gj=5.0/1000;#gene concentration uM
  RXT=0.15;#RNAP concentration uM
  RLT=1.6;#ribosome concentration uM
  KX=0.3;#transcription saturation constant uM
  KL=57.0;#translation saturation constant uM
  tauX=2.7;#transcription time constant
  tauL=0.8;#translation time constant
  kdX=8.35/3600;#mRNA degradation rate s^-1
  kdL=(9.9e-3)/3600;#protein degradation rate s^-1

#We have to solve the FBA at each different inducer concentration
  I=collect(0.0001:0.0001:10.0);#mM
  translation_rate=zeros(length(I));
#Run FBA for every inducer concentration
for i in 1:length(I)
  #Create transcription rate equation
  rX = kEX*RXT*Gj/(KX*tauX+(tauX+1)*Gj);
  #Create control function dependent on weights and binding paramters
  W1=0.26;
  W2=300.0;
  K=0.3;#mM
  n=1.5;
  fI = (I[i]^n)/((I[i]^n)+K^n);
  u = (W1+W2*fI)/(1+W1+W2*fI);
  r_hatX = rX*u;
  #Create translation rate equation
  mRNA_steady_state = r_hatX/kdX;#the mRNA value is at steady state since this value is set constant
  r_hatL=kEL*RLT*mRNA_steady_state/(KL*tauL+(tauL+1)*mRNA_steady_state);

#Create flux bounds
global Default_bounds_array = zeros(15,2);
#Set upper bounds
Default_bounds_array[1,2]=Inf; #v1 transcription initiation. Not given
Default_bounds_array[2,2]=r_hatX; #v2 transcription
Default_bounds_array[3,2]=mRNA_steady_state*kdX; #v3 mRNA decay
Default_bounds_array[4,2]=Inf; #v4 translation initiation. Not given
Default_bounds_array[5,2]=r_hatL; #v5 translation
Default_bounds_array[6,2]=Inf; #v6 tRNA charging. Not given
Default_bounds_array[7,2]=100000/3600; #b1 AA influx um/s
Default_bounds_array[8,2]=100000/3600; #b2 NTP influx
Default_bounds_array[9,2]=100000/3600; #b3 protein outflux
Default_bounds_array[10,2]=100000/3600; #b4 NMP outflux
Default_bounds_array[11,2]=100000/3600; #b5 ATP influx
Default_bounds_array[12,2]=100000/3600; #b6 AMP outflux
Default_bounds_array[13,2]=100000/3600; #b7 GTP influx
Default_bounds_array[14,2]=100000/3600; #b8 GDP outflux
Default_bounds_array[15,2]=100000/3600; #b9 Pi outflux
#Set lower bounds
Default_bounds_array[1,1]=0.0; #v1 transcription initiation. Nonreversible
Default_bounds_array[2,1]=r_hatX; #v2 transcription (value is constant)
Default_bounds_array[3,1]=mRNA_steady_state*kdX; #v3 mRNA decay (value is constant)
Default_bounds_array[4,1]=0.0; #v4 translation initiation. Nonreversible
Default_bounds_array[5,1]=0.0; #v5 translation. Nonreversible
Default_bounds_array[6,1]=0.0; #v6 tRNA charging. Nonreversible
Default_bounds_array[7,1]=-100000/3600; #b1 AA influx um/s
Default_bounds_array[8,1]=-100000/3600; #b2 NTP influx
Default_bounds_array[9,1]=-100000/3600; #b3 protein outflux
Default_bounds_array[10,1]=-100000/3600; #b4 NMP outflux
Default_bounds_array[11,1]=-100000/3600; #b5 ATP influx
Default_bounds_array[12,1]=-100000/3600; #b6 AMP outflux
Default_bounds_array[13,1]=-100000/3600; #b7 GTP influx
Default_bounds_array[14,1]=-100000/3600; #b8 GDP outflux
Default_bounds_array[15,1]=-100000/3600; #b9 Pi outflux

#Create species bounds arrays (no bounds)
Species_bounds_array=zeros(17,2);

#Create objective function coefficient array. Only translation (v5) is maximized:
Objective_coefficient_array=zeros(15);
Objective_coefficient_array[5]=-1;

Min_flag = true;

#Use the FBA function
answer=calculate_optimal_flux_distribution(Stoichiometric_matrix,Default_bounds_array,Species_bounds_array,Objective_coefficient_array);
r=answer[2];
translation_rate[i]=r[5];
end

#Extract the steady state protein level from the translation rate by dividing
  #by the degradation rate (same method as for mRNA)
protein_steady_state=translation_rate./kdL;
I_plot=log10.(I);
plot(I_plot,protein_steady_state)
title("Protein as a Function of Inducer")
ylabel("Protein uM")
xlabel("Inducer log(mM)")
