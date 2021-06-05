using CSV

Busdata=CSV.read("Data/Busdata.csv", DataFrame);
Gendata=CSV.read("Data/Gendata.csv", DataFrame);
Gencostdata=CSV.read("Data/Gencostdata.csv", DataFrame);
Branchdata=CSV.read("Data/Branchdata.csv", DataFrame);
Linedata_candidate=CSV.read("Data/Linedata_candidate.csv", DataFrame);


#bus Pd Qd 
Loaddata=zeros(1,3);
for i in 1:size(Busdata,1)
    if Busdata[i,3] !=0
        global Loaddata=[Loaddata;i 1*Busdata[i,3] 1*Busdata[i,4]]; #Load level modification
    end
end
Loaddata=Loaddata[2:end,:];

Num_bus=size(Busdata,1);
Num_gen=size(Gendata,1);
Num_load=count(!iszero, Busdata[:,3]);
Num_Eline=size(Branchdata,1);
Num_Cline=size(Linedata_candidate,1);
Loadshed_cost=9000#($/MWh);
Maxinvest_cost=2000000#(M$);
baseMVA=100;
ELcap=1; #Exisitng lines capacity reduction factor

ELcap=1; #Exisitng lines capacity reduction factor
GE_ff=zeros(Num_Eline);GE_ft=zeros(Num_Eline);GE_tf=zeros(Num_Eline);GE_tt=zeros(Num_Eline);
BE_ff=zeros(Num_Eline);BE_ft=zeros(Num_Eline);BE_tf=zeros(Num_Eline);BE_tt=zeros(Num_Eline);
for l=1:Num_Eline
    if Branchdata[l,9]==0
       GE_ff[l]=Branchdata[l,3]/(Branchdata[l,3]^2+Branchdata[l,4]^2)
       GE_ft[l]=-Branchdata[l,3]/(Branchdata[l,3]^2+Branchdata[l,4]^2)
       GE_tf[l]=GE_ft[l]
       GE_tt[l]=GE_ff[l]
       BE_ff[l]=(Branchdata[l,5]/2)-Branchdata[l,4]/(Branchdata[l,3]^2+Branchdata[l,4]^2)
       BE_ft[l]=Branchdata[l,4]/(Branchdata[l,3]^2+Branchdata[l,4]^2)
       BE_tf[l]=BE_ft[l]
       BE_tt[l]=BE_ff[l]
    else
       GE_ff[l]=(Branchdata[l,9]^2)*(Branchdata[l,3]/(Branchdata[l,3]^2+Branchdata[l,4]^2))
       GE_ft[l]=(-Branchdata[l,9])*(Branchdata[l,3]/(Branchdata[l,3]^2+Branchdata[l,4]^2))
       GE_tf[l]=GE_ft[l]
       GE_tt[l]=Branchdata[l,3]/(Branchdata[l,3]^2+Branchdata[l,4]^2)
       BE_ff[l]=(Branchdata[l,9]^2)*(-Branchdata[l,4]/(Branchdata[l,3]^2+Branchdata[l,4]^2))
       BE_ft[l]=(-Branchdata[l,9])*(-Branchdata[l,4]/(Branchdata[l,3]^2+Branchdata[l,4]^2))
       BE_tf[l]=BE_ft[l]
       BE_tt[l]=(-Branchdata[l,4]/(Branchdata[l,3]^2+Branchdata[l,4]^2))
    end
end

GC_ff=zeros(Num_Cline);GC_ft=zeros(Num_Cline);GC_tf=zeros(Num_Cline);GC_tt=zeros(Num_Cline);
BC_ff=zeros(Num_Cline);BC_ft=zeros(Num_Cline);BC_tf=zeros(Num_Cline);BC_tt=zeros(Num_Cline);
for l=1:Num_Cline
    GC_ff[l]=Linedata_candidate[l,3]/(Linedata_candidate[l,3]^2+Linedata_candidate[l,4]^2)
    GC_ft[l]=-Linedata_candidate[l,3]/(Linedata_candidate[l,3]^2+Linedata_candidate[l,4]^2)
    GC_tf[l]=GC_ft[l]
    GC_tt[l]=GC_ff[l]
    BC_ff[l]=-Linedata_candidate[l,4]/(Linedata_candidate[l,3]^2+Linedata_candidate[l,4]^2)
    BC_ft[l]=Linedata_candidate[l,4]/(Linedata_candidate[l,3]^2+Linedata_candidate[l,4]^2)
    BC_tf[l]=BC_ft[l]
    BC_tt[l]=BC_ff[l]
end
     
Gsh=zeros(Num_bus);Bsh=zeros(Num_bus);
Gsh=Busdata[:,5];Bsh=Busdata[:,6];

Theta_min=-(pi/2)*ones(Num_Eline+Num_Cline);
Theta_max=(pi/2)*ones(Num_Eline+Num_Cline);
RE_min=zeros(Num_Eline);
RE_max=1.05*1.05*ones(Num_Eline);
LE_min=-1.05*1.05*ones(Num_Eline);
LE_max=1.05*1.05*ones(Num_Eline);
V_min=0.95*ones(Num_bus);
V_max=1.05*ones(Num_bus);
CR_min=zeros(Num_Eline);
CR_max=ones(Num_Eline);
SR_min=-ones(Num_Eline);
SR_max=ones(Num_Eline);


struct Bounds
    Theta_min
    V_min
    RE_min
    LE_min
    CR_min
    SR_min
    Theta_max
    V_max
    RE_max
    LE_max
    CR_max
    SR_max
end

bounds = Bounds(Theta_min, V_min, RE_min, LE_min, CR_min, SR_min, Theta_max, V_max, RE_max, LE_max, CR_max, SR_max);

XL=zeros(Num_Cline)


