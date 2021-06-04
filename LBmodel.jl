function relax(bounds,RT,K)

Theta_min = bounds.Theta_min
V_min = bounds.V_min
RE_min = bounds.RE_min
LE_min = bounds.LE_min
CR_min = bounds.CR_min
SR_min = bounds.SR_min
Theta_max = bounds.Theta_max
V_max = bounds.V_max
RE_max = bounds.RE_max
LE_max = bounds.LE_max
CR_max = bounds.CR_max
SR_max = bounds.SR_max

model=Model(CPLEX.Optimizer)
set_optimizer_attribute(model, "CPX_PARAM_EPGAP", 0.001)#Relative MIP gap tolerance

@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=PLe_S[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=QLe_S[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=PLe_R[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=QLe_R[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=PLc_S[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=QLc_S[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=PLc_R[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=QLc_R[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,0<=Pg[i in 1:Num_gen]<=Gendata[i,9], start=Gendata[i,9])
@variable(model,Gendata[i,5]<=Qg[i in 1:Num_gen]<=Gendata[i,4], start=Gendata[i,4])
@variable(model,XL[1:Num_Cline],Bin)
@variable(model,0<=Qi[1:Num_bus]<50)

#original variables
@variable(model,V_min[i]<=Vi[i in 1:Num_bus]<=V_max[i], start=1)
@variable(model,-pi/2<=theta[1:Num_bus]<=pi/2, start=0)
@variable(model,Theta_min[i]<=Theta[i in 1:Num_Eline]<=Theta_max[i], start=0)

#auxiliary variables
@variable(model,V_min[i]^2<=Ui[i in 1:Num_bus]<=V_max[i]^2, start=1)
@variable(model,RE_min[i]<=RE[i in 1:Num_Eline]<=RE_max[i], start=1)
@variable(model,-1.05*1.05<=RC[i in 1:Num_Cline]<=1.05*1.05, start=1)
@variable(model,LE_min[i]<=LE[i in 1:Num_Eline]<=LE_max[i], start=0)
@variable(model,-1.05*1.05<=LC[i in 1:Num_Cline]<=1.05*1.05, start=0)

@constraint(model,con_budgetconstraint,sum(Linedata_candidate[i,6]*XL[i] for i in 1:Num_Cline)<=Maxinvest_cost)
@constraint(model,con_refbus,sum(theta[i] for i in 1:Num_bus if Busdata[i,2]==3)==0)
@constraint(model,con_Thetaexisitng[i=1:Num_Eline],Theta[i]-(theta[Branchdata[i,1]]-theta[Branchdata[i,2]])==0)

@constraint(model,con_nodalActivepowerbalance[i=1:Num_bus],-sum(Pg[j] for j in 1:Num_gen if Gendata[j,1]==i)+
baseMVA*sum(PLe_R[j] for j in 1:Num_Eline if Branchdata[j,2]==i)+baseMVA*sum(PLc_R[j] for j in 1:Num_Cline if Linedata_candidate[j,2]==i)+
baseMVA*sum(PLe_S[j] for j in 1:Num_Eline if Branchdata[j,1]==i)+baseMVA*sum(PLc_S[j] for j in 1:Num_Cline if Linedata_candidate[j,1]==i)+
sum(Loaddata[j,2] for j in 1:Num_load if Loaddata[j,1]==i)+Gsh[i]*Ui[i]==0)
@constraint(model,con_nodalReactivepowerbalance[i=1:Num_bus],-sum(Qg[j] for j in 1:Num_gen if Gendata[j,1]==i)+
baseMVA*sum(QLe_R[j] for j in 1:Num_Eline if Branchdata[j,2]==i)+baseMVA*sum(QLc_R[j] for j in 1:Num_Cline if Linedata_candidate[j,2]==i)+
baseMVA*sum(QLe_S[j] for j in 1:Num_Eline if Branchdata[j,1]==i)+baseMVA*sum(QLc_S[j] for j in 1:Num_Cline if Linedata_candidate[j,1]==i)+
sum(Loaddata[j,3] for j in 1:Num_load if Loaddata[j,1]==i)-Bsh[i]*Ui[i]-Qi[i]==0)

#Exisitng line constraints
@constraint(model,con_ExisitngActivelineflowSend[i=1:Num_Eline],-PLe_S[i]+GE_ff[i]*Ui[Branchdata[i,1]]+GE_ft[i]*RE[i]+BE_ft[i]*LE[i]==0)
@constraint(model,con_ExisitngReactivelineflowSenc[i=1:Num_Eline],-QLe_S[i]-BE_ff[i]*Ui[Branchdata[i,1]]-BE_ft[i]*RE[i]+GE_ft[i]*LE[i]==0)
@constraint(model,con_ExisitngActivelineflowRec[i=1:Num_Eline],-PLe_R[i]+GE_tt[i]*Ui[Branchdata[i,2]]+GE_tf[i]*RE[i]-BE_tf[i]*LE[i]==0)
@constraint(model,con_ExisitngReactivelineflowRec[i=1:Num_Eline],-QLe_R[i]-BE_tt[i]*Ui[Branchdata[i,2]]-BE_tf[i]*RE[i]-GE_tf[i]*LE[i]==0)

#Candidate line constraints
@constraint(model,con_CandidateActivelineflowmaxSend[i=1:Num_Cline],-PLc_S[i]+GC_ff[i]*Ui[convert(UInt8,Linedata_candidate[i,1])]+GC_ft[i]*RC[i]+BC_ft[i]*LC[i]-
((abs(GC_ff[i])+abs(GC_ft[i])+abs(BC_ft[i]))*1.05^2)*(1-XL[i])<=0)
@constraint(model,con_CandidateActivelineflowminSend[i=1:Num_Cline],-PLc_S[i]+GC_ff[i]*Ui[convert(UInt8,Linedata_candidate[i,1])]+GC_ft[i]*RC[i]+BC_ft[i]*LC[i]+
((abs(GC_ff[i])+abs(GC_ft[i])+abs(BC_ft[i]))*1.05^2)*(1-XL[i])>=0)
@constraint(model,con_LineActiveflowlmax_candidateSend[i=1:Num_Cline],PLc_S[i]-XL[i]*(Linedata_candidate[i,5]/baseMVA)<=0)
@constraint(model,con_LineActiveflowmin_candidateSend[i=1:Num_Cline],PLc_S[i]+XL[i]*(Linedata_candidate[i,5]/baseMVA)>=0)

@constraint(model,con_CandidateReactivelineflowmaxSend[i=1:Num_Cline],-QLc_S[i]-BC_ff[i]*Ui[convert(UInt8,Linedata_candidate[i,1])]-BC_ft[i]*RC[i]+GC_ft[i]*LC[i]-
((abs(BC_ff[i])+abs(BC_ft[i])+abs(GC_ft[i]))*1.05^2)*(1-XL[i])<=0)
@constraint(model,con_CandidateReactivelineflowminSend[i=1:Num_Cline],-QLc_S[i]-BC_ff[i]*Ui[convert(UInt8,Linedata_candidate[i,1])]-BC_ft[i]*RC[i]+GC_ft[i]*LC[i]+
((abs(BC_ff[i])+abs(BC_ft[i])+abs(GC_ft[i]))*1.05^2)*(1-XL[i])>=0)
@constraint(model,con_LineReactiveflowlmax_candidateSend[i=1:Num_Cline],QLc_S[i]-XL[i]*(Linedata_candidate[i,5]/baseMVA)<=0)
@constraint(model,con_LineReactiveflowmin_candidateSend[i=1:Num_Cline],QLc_S[i]+XL[i]*(Linedata_candidate[i,5]/baseMVA)>=0)

@constraint(model,con_CandidateActivelineflowmaxRec[i=1:Num_Cline],-PLc_R[i]+GC_tt[i]*Ui[convert(UInt8,Linedata_candidate[i,2])]+GC_tf[i]*RC[i]-BC_tf[i]*LC[i]-
((abs(GC_tt[i])+abs(GC_tf[i])+abs(BC_tf[i]))*1.05^2)*(1-XL[i])<=0)
@constraint(model,con_CandidateActivelineflowminRec[i=1:Num_Cline],-PLc_R[i]+GC_tt[i]*Ui[convert(UInt8,Linedata_candidate[i,2])]+GC_tf[i]*RC[i]-BC_tf[i]*LC[i]+
((abs(GC_tt[i])+abs(GC_tf[i])+abs(BC_tf[i]))*1.05^2)*(1-XL[i])>=0)
@constraint(model,con_LineActiveflowlmax_candidateRec[i=1:Num_Cline],PLc_R[i]-XL[i]*(Linedata_candidate[i,5]/baseMVA)<=0)
@constraint(model,con_LineActiveflowmin_candidateRec[i=1:Num_Cline],PLc_R[i]+XL[i]*(Linedata_candidate[i,5]/baseMVA)>=0)

@constraint(model,con_CandidateReactivelineflowmaxRec[i=1:Num_Cline],-QLc_R[i]-BC_tt[i]*Ui[convert(UInt8,Linedata_candidate[i,2])]-BC_tf[i]*RC[i]-GC_tf[i]*LC[i]-
((abs(BC_tt[i])+abs(BC_tf[i])+abs(GC_tf[i]))*1.05^2)*(1-XL[i])<=0)
@constraint(model,con_CandidateReactivelineflowminRec[i=1:Num_Cline],-QLc_R[i]-BC_tt[i]*Ui[convert(UInt8,Linedata_candidate[i,2])]-BC_tf[i]*RC[i]-GC_tf[i]*LC[i]+
((abs(BC_tt[i])+abs(BC_tf[i])+abs(GC_tf[i]))*1.05^2)*(1-XL[i])>=0)
@constraint(model,con_LineReactiveflowlmax_candidateRec[i=1:Num_Cline],QLc_R[i]-XL[i]*(Linedata_candidate[i,5]/baseMVA)<=0)
@constraint(model,con_LineReactiveflowmin_candidateRec[i=1:Num_Cline],QLc_R[i]+XL[i]*(Linedata_candidate[i,5]/baseMVA)>=0)

#Thermal constraint--------------------------------------------------------------------------
@constraint(model,con_thermalElinelimitSend[i=1:Num_Eline],(PLe_S[i])^2+(QLe_S[i])^2<=(ELcap*Branchdata[i,6]/baseMVA)^2)
@constraint(model,con_thermalElinelimitRec[i=1:Num_Eline],(PLe_R[i])^2+(QLe_R[i])^2<=(ELcap*Branchdata[i,6]/baseMVA)^2)

#Conic constraint--------------------------------------------------------------------------
@constraint(model,con_conic1[i=1:Num_Eline],RE[i]^2+LE[i]^2<=Ui[Branchdata[i,1]]*Ui[Branchdata[i,2]])
@constraint(model,con_conic2[i=1:Num_Cline],RC[i]^2+LC[i]^2<=Ui[convert(UInt8,Linedata_candidate[i,1])]*Ui[convert(UInt8,Linedata_candidate[i,2])])

if RT>0
   
   #Relaxation tightening-----------------------------------------------------------------------------------------------------------
   Flag_thermal=ones(Num_Eline);Flag_thermal[1]=1;
   Flag_cycle=ones(100000);Flag_cycle[1]=1;
   Flag_conic=ones(Num_Eline);Flag_conic[1]=1;
   Flag_quad=ones(Num_bus);Flag_quad[1]=1;
   
      
   #Cycle constraint--------------------------------------------------------------------------
   #atan constraint------
   SLU=zeros(Num_Eline);SUU=zeros(Num_Eline);SLO=zeros(Num_Eline);SUO=zeros(Num_Eline);
   for i=1:Num_Eline
       function f!(F, x)
            F[1]=((atan(-LE_min[i]/RE_min[i])-atan(x[1]/RE_min[i]))/(-LE_min[i]-x[1]))-(RE_min[i]/(x[1]^2+RE_min[i]^2))
       end
       rrr=nlsolve(f!,[-LE_max[i]])
       SLU[i]=max(-LE_max[i],rrr.zero[1])
       
       function g!(F, x)
            F[1]=((atan(x[1]/RE_max[i])-atan(-LE_max[i]/RE_min[i]))/(x[1]+LE_max[i]))-(RE_max[i]/(LE_max[i]^2+RE_max[i]^2))
       end
       rrr=nlsolve(g!,[-LE_min[i]])
       SUU[i]=max(-LE_min[i],rrr.zero[1])
   
       function h!(F, x)
            F[1]=((atan(-LE_min[i]/RE_max[i])-atan(x[1]/RE_max[i]))/(-LE_min[i]-x[1]))-(RE_max[i]/(LE_min[i]^2+RE_max[i]^2))
       end
       rrr=nlsolve(h!,[-LE_max[i]])
       SLO[i]=min(-LE_max[i],rrr.zero[1])
   
       function I!(F, x)
            F[1]=((atan(x[1]/RE_min[i])-atan(-LE_max[i]/RE_min[i]))/(x[1]+LE_max[i]))-(RE_min[i]/(x[1]^2+RE_min[i]^2))
       end
       rrr=nlsolve(I!,[-LE_min[i]])
       SUO[i]=min(-LE_min[i],rrr.zero[1])
   end
   
   alpha_u=zeros(Num_Eline,2);alpha_o=zeros(Num_Eline,2);beta_u=zeros(Num_Eline,2);beta_o=zeros(Num_Eline,2);
   for i=1:Num_Eline
       alpha_u[i,1]=(atan(-LE_min[i]/RE_min[i])-atan(SLU[i]/RE_min[i]))/(-LE_min[i]-SLU[i])
       alpha_u[i,2]=(atan(SUU[i]/RE_max[i])-atan(-LE_max[i]/RE_max[i]))/(SUU[i]+LE_max[i])
       alpha_o[i,1]=(atan(-LE_min[i]/RE_max[i])-atan(SLO[i]/RE_max[i]))/(-LE_min[i]-SLO[i])
       alpha_o[i,2]=(atan(SUO[i]/RE_min[i])-atan(-LE_max[i]/RE_min[i]))/(SUO[i]+LE_max[i])
       beta_u[i,1]=(LE_min[i])/(LE_min[i]^2+RE_min[i]^2)
       beta_u[i,2]=(atan(-LE_max[i]/RE_max[i])-atan(-LE_max[i]/RE_min[i]))/(RE_max[i]-RE_min[i])
       beta_o[i,1]=(atan(-LE_min[i]/RE_max[i])-atan(-LE_min[i]/RE_min[i]))/(RE_max[i]-RE_min[i])
       beta_o[i,2]=(LE_max[i])/(LE_max[i]^2+RE_min[i]^2)
   end
   
   gamma_u=zeros(Num_Eline,2);gamma_o=zeros(Num_Eline,2);
   for i=1:Num_Eline
       gamma_u[i,1]=max((alpha_u[i,1]*(-LE_min[i])+beta_u[i,1]*(RE_min[i]))-atan(-LE_min[i]/RE_min[i]),(alpha_u[i,1]*(-LE_max[i])+beta_u[i,1]*(RE_min[i]))-atan(-LE_max[i]/RE_min[i]),
                        (alpha_u[i,1]*(-LE_min[i])+beta_u[i,1]*(RE_max[i]))-atan(-LE_min[i]/RE_max[i]),(alpha_u[i,1]*(-LE_max[i])+beta_u[i,1]*(RE_max[i]))-atan(-LE_max[i]/RE_max[i]))
       gamma_u[i,2]=max((alpha_u[i,2]*(-LE_min[i])+beta_u[i,2]*(RE_min[i]))-atan(-LE_min[i]/RE_min[i]),(alpha_u[i,2]*(-LE_max[i])+beta_u[i,2]*(RE_min[i]))-atan(-LE_max[i]/RE_min[i]),
                        (alpha_u[i,2]*(-LE_min[i])+beta_u[i,2]*(RE_max[i]))-atan(-LE_min[i]/RE_max[i]),(alpha_u[i,2]*(-LE_max[i])+beta_u[i,2]*(RE_max[i]))-atan(-LE_max[i]/RE_max[i]))
       gamma_o[i,1]=max(atan(-LE_min[i]/RE_min[i])-(alpha_o[i,1]*(-LE_min[i])+beta_o[i,1]*(RE_min[i])),atan(-LE_max[i]/RE_min[i])-(alpha_o[i,1]*(-LE_max[i])+beta_o[i,1]*(RE_min[i])),
                        atan(-LE_min[i]/RE_max[i])-(alpha_o[i,1]*(-LE_min[i])+beta_o[i,1]*(RE_max[i])),atan(-LE_max[i]/RE_max[i])-(alpha_o[i,1]*(-LE_max[i])+beta_o[i,1]*(RE_max[i])))
       gamma_o[i,2]=max(atan(-LE_min[i]/RE_min[i])-(alpha_o[i,2]*(-LE_min[i])+beta_o[i,2]*(RE_min[i])),atan(-LE_max[i]/RE_min[i])-(alpha_o[i,2]*(-LE_max[i])+beta_o[i,2]*(RE_min[i])),
                        atan(-LE_min[i]/RE_max[i])-(alpha_o[i,2]*(-LE_min[i])+beta_o[i,2]*(RE_max[i])),atan(-LE_max[i]/RE_max[i])-(alpha_o[i,2]*(-LE_max[i])+beta_o[i,2]*(RE_max[i])))
   end
   
   @constraint(model,[i=1:Num_Eline,k=1:2],Theta[i]+alpha_u[i,k]*(-LE[i])+beta_u[i,k]*(RE[i])-gamma_u[i,k]<=0)
   @constraint(model,[i=1:Num_Eline,k=1:2],Theta[i]+alpha_o[i,k]*(-LE[i])+beta_o[i,k]*(RE[i])+gamma_o[i,k]>=0)
   
   
   #--------------------
   @constraint(model,con_cyc1[i=1:Num_Eline],LE[i]-RE[i]*tan(Theta_max[i])<=0)
   @constraint(model,con_cyc2[i=1:Num_Eline],LE[i]-RE[i]*tan(Theta_min[i])>=0)
   g=SimpleDiGraph(Num_bus)
   for l=1:Num_Eline
       add_edge!(g, Branchdata[l,1], Branchdata[l,2])
   end
   g=union(g,reverse(g))
   C=filter(x->length(x)>=3,cycle_basis(g))
   for i=1:size(C,1)
       cyc=C[i]
       edge=zeros(size(cyc,1))
       for j=1:size(cyc,1)-1
           for l=1:Num_Eline
               if (Branchdata[l,1]==cyc[j]) & (Branchdata[l,2]==cyc[j+1])
                  edge[j]=l
               elseif (Branchdata[l,2]==cyc[j]) & (Branchdata[l,1]==cyc[j+1])
                  edge[j]=-l
               end
           end 
        end
        for l=1:Num_Eline
            if (Branchdata[l,1]==cyc[end]) & (Branchdata[l,2]==cyc[1])
               edge[end]=l
            elseif (Branchdata[l,2]==cyc[end]) & (Branchdata[l,1]==cyc[1])
               edge[end]=-l
            end
        end
        if Flag_cycle[i] !=0        
           @constraint(model,sum(sign(edge[h])*Theta[convert(UInt8,edge[h]*sign(edge[h]))] for h=1:size(cyc,1) )==0)
           for h=1:size(cyc,1) 
               Flag_conic[convert(UInt8,edge[h]*sign(edge[h]))]=1
           end
        end
   end
   
   #Cosine and Sine constraints-------------------------------------------------------------------------
   
   for l=1:Num_Eline
       if (Theta_max[l]>0) & (Theta_min[l]>0)
          CR_min[l]=max(CR_min[l],min(cos(Theta_max[l]),cos(Theta_min[l])));
          CR_max[l]=min(CR_max[l],max(cos(Theta_max[l]),cos(Theta_min[l])));
          SR_min[l]=max(SR_min[l],min(sin(Theta_max[l]),sin(Theta_min[l])));
          SR_max[l]=min(SR_max[l],max(sin(Theta_max[l]),sin(Theta_min[l])));
       elseif (Theta_max[l]<0) & (Theta_min[l]<0)
          CR_min[l]=max(CR_min[l],min(cos(Theta_max[l]),cos(Theta_min[l])));
          CR_max[l]=min(CR_max[l],max(cos(Theta_max[l]),cos(Theta_min[l])));
          SR_min[l]=max(SR_min[l],min(sin(Theta_max[l]),sin(Theta_min[l])));
          SR_max[l]=min(SR_max[l],max(sin(Theta_max[l]),sin(Theta_min[l])));
       elseif (Theta_max[l]>0) & (Theta_min[l]<0)
          CR_min[l]=max(CR_min[l],min(cos(Theta_max[l]),cos(Theta_min[l])));
          CR_max[l]=min(CR_max[l],1);
          SR_min[l]=max(SR_min[l],sin(Theta_min[l]));
          SR_max[l]=min(SR_max[l],sin(Theta_max[l]));
       end
   end
   
   @variable(model,CR_min[i]<=CR[i in 1:Num_Eline]<=CR_max[i], start=0)
   @variable(model,SR_min[i]<=SR[i in 1:Num_Eline]<=SR_max[i], start=0)
   
   Theta_M=zeros(Num_Eline);
   Theta_u=10*ones(Num_Eline);
   r=zeros(Num_Eline);
   for l=1:Num_Eline
       Theta_M[l]=max(abs(Theta_max[l]),abs(Theta_min[l]));
       Theta_u[l]=(Theta_max[l]+Theta_min[l])/2;
       if (Theta_max[l]>0) & (Theta_min[l]>0)
          r[l]=1;
       elseif (Theta_max[l]<0) & (Theta_min[l]<0)
          r[l]=-1;
       elseif (Theta_max[l]>0) & (Theta_min[l]<0)
          r[l]=0;
       end
   end
   
   
   @constraint(model,con1_CR[i=1:Num_Eline;Flag_conic[i] != 0],CR[i]-1+((1-cos(Theta_M[i]))/Theta_M[i]^2)*Theta[i]^2<=0)
   @constraint(model,con2_CR[i=1:Num_Eline;Flag_conic[i] != 0],CR[i]-((cos(Theta_max[i])-cos(Theta_min[i]))/(Theta_max[i]-Theta_min[i]))*(Theta[i]-Theta_max[i])-cos(Theta_max[i])>=0)
   
   for i=1:Num_Eline
       if r[i]==0 & (Flag_conic[i] != 0)
           @constraint(model,SR[i]-cos(Theta_M[i]/2)*(Theta[i]-Theta_M[i]/2)-sin(Theta_M[i]/2)<=0)
           @constraint(model,SR[i]-cos(Theta_M[i]/2)*(Theta[i]+Theta_M[i]/2)+sin(Theta_M[i]/2)>=0)
       elseif r[i]==1 & (Flag_conic[i] != 0)
       	@constraint(model,SR[i]-((sin(Theta_max[i])-sin(Theta_min[i]))/(Theta_max[i]-Theta_min[i]))*(Theta[i]-Theta_min[i])-sin(Theta_min[i])>=0)
   	@constraint(model,SR[i]-sin(Theta_min[i])-cos(Theta_min[i])*(Theta[i]-Theta_min[i])<=0)
   	@constraint(model,SR[i]-sin(Theta_max[i])-cos(Theta_max[i])*(Theta[i]-Theta_max[i])<=0)
   	@constraint(model,SR[i]-sin(Theta_u[i])-cos(Theta_u[i])*(Theta[i]-Theta_u[i])<=0)
       elseif r[i]==-1 & (Flag_conic[i] != 0)
       	@constraint(model,SR[i]-((sin(Theta_max[i])-sin(Theta_min[i]))/(Theta_max[i]-Theta_min[i]))*(Theta[i]-Theta_min[i])-sin(Theta_min[i])<=0)
   	@constraint(model,SR[i]-sin(Theta_min[i])-cos(Theta_min[i])*(Theta[i]-Theta_min[i])>=0)
   	@constraint(model,SR[i]-sin(Theta_max[i])-cos(Theta_max[i])*(Theta[i]-Theta_max[i])>=0)
   	@constraint(model,SR[i]-sin(Theta_u[i])-cos(Theta_u[i])*(Theta[i]-Theta_u[i])>=0)
       end
   end   
   
   if K==1
      
      Num_pieces=2
      
      #Piecewise quadratic constraints as lambda------------------------------------------------------------------
      KK=Num_pieces+1;
      Xi_V=zeros(Num_bus,KK);
      for i=1:Num_bus
          for m=1:Num_pieces
              Xi_V[i,m]=V_min[i]+((m-1)/Num_pieces).*(V_max[i]-V_min[i]);
          end
      end
      Xi_V[:,KK]=V_max[:];
      
      @variable(model,0<=lq[1:Num_bus,1:KK])
      @variable(model,zq[1:Num_bus,1:Num_pieces],Bin)
      
      @constraint(model,[i=1:Num_bus],sum(zq[i,j] for j=1:Num_pieces)==1)
      @constraint(model,[i=1:Num_bus],sum(lq[i,j] for j=1:KK)==1)
      
      @constraint(model,[i=1:Num_bus],zq[i,1]-lq[i,1]>=0)
      @constraint(model,[i=1:Num_bus],zq[i,KK-1]-lq[i,KK]>=0)
      @constraint(model,[i=1:Num_bus,k=2:KK-1],zq[i,k]+zq[i,k-1]-lq[i,k]>=0)
      
      @constraint(model,[i=1:Num_bus],Vi[i]-sum(lq[i,k]*Xi_V[i,k] for k=1:KK)==0)
      @constraint(model,[i=1:Num_bus],Ui[i]-Vi[i]^2>=0)
      @constraint(model,[i=1:Num_bus],Ui[i]-sum(lq[i,k]*(Xi_V[i,k])^2 for k=1:KK)<=0)
      
      #Strengthening valid inequalities-------------------------------------
      @constraint(model,[i=1:Num_bus],Vi[i]-sum(zq[i,k-1]*Xi_V[i,k] for k=2:KK)<=0)
      @constraint(model,[i=1:Num_bus],Vi[i]-sum(zq[i,k]*Xi_V[i,k] for k=1:KK-1)>=0)
      
           
      #Trilinear constraints as lambda--------------------------------------{L=Vi*Vj*SRij}------------------------------
      @variable(model,0<=lR[1:Num_Eline,1:8])
      @variable(model,0<=lL[1:Num_Eline,1:8])
      
      @constraint(model,[i=1:Num_Eline],sum(lR[i,j] for j=1:8)==1)
      @constraint(model,[i=1:Num_Eline],sum(lL[i,j] for j=1:8)==1)
      
      @constraint(model,[i=1:Num_Eline],RE[i]-
      (lR[i,1]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR_min[i]+
      lR[i,2]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR_max[i]+
      lR[i,3]*V_min[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR_min[i]+
      lR[i,4]*V_min[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR_max[i]+
      lR[i,5]*V_max[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR_min[i]+
      lR[i,6]*V_max[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR_max[i]+
      lR[i,7]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR_min[i]+
      lR[i,8]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR_max[i])==0)
      
      @constraint(model,[i=1:Num_Eline],Vi[Branchdata[i,1]]-(lR[i,1]*V_min[Branchdata[i,1]]+lR[i,2]*V_min[Branchdata[i,1]]+lR[i,3]*V_min[Branchdata[i,1]]+lR[i,4]*V_min[Branchdata[i,1]]+
      lR[i,5]*V_max[Branchdata[i,1]]+lR[i,6]*V_max[Branchdata[i,1]]+lR[i,7]*V_max[Branchdata[i,1]]+lR[i,8]*V_max[Branchdata[i,1]])==0)
      @constraint(model,[i=1:Num_Eline],Vi[Branchdata[i,2]]-(lR[i,1]*V_min[Branchdata[i,2]]+lR[i,2]*V_min[Branchdata[i,2]]+lR[i,3]*V_max[Branchdata[i,2]]+lR[i,4]*V_max[Branchdata[i,2]]+
      lR[i,5]*V_min[Branchdata[i,2]]+lR[i,6]*V_min[Branchdata[i,2]]+lR[i,7]*V_max[Branchdata[i,2]]+lR[i,8]*V_max[Branchdata[i,2]])==0)
      @constraint(model,[i=1:Num_Eline],CR[i]-(lR[i,1]*CR_min[i]+lR[i,2]*CR_max[i]+lR[i,3]*CR_min[i]+lR[i,4]*CR_max[i]+
      lR[i,5]*CR_min[i]+lR[i,6]*CR_max[i]+lR[i,7]*CR_min[i]+lR[i,8]*CR_max[i])==0)
      
      
      @constraint(model,[i=1:Num_Eline],LE[i]-
      (lL[i,1]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR_min[i]+
      lL[i,2]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR_max[i]+
      lL[i,3]*V_min[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR_min[i]+
      lL[i,4]*V_min[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR_max[i]+
      lL[i,5]*V_max[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR_min[i]+
      lL[i,6]*V_max[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR_max[i]+
      lL[i,7]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR_min[i]+
      lL[i,8]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR_max[i])==0)
      
      @constraint(model,[i=1:Num_Eline],Vi[Branchdata[i,1]]-(lL[i,1]*V_min[Branchdata[i,1]]+lL[i,2]*V_min[Branchdata[i,1]]+lL[i,3]*V_min[Branchdata[i,1]]+lL[i,4]*V_min[Branchdata[i,1]]+
      lL[i,5]*V_max[Branchdata[i,1]]+lL[i,6]*V_max[Branchdata[i,1]]+lL[i,7]*V_max[Branchdata[i,1]]+lL[i,8]*V_max[Branchdata[i,1]])==0)
      @constraint(model,[i=1:Num_Eline],Vi[Branchdata[i,2]]-(lL[i,1]*V_min[Branchdata[i,2]]+lL[i,2]*V_min[Branchdata[i,2]]+lL[i,3]*V_max[Branchdata[i,2]]+lL[i,4]*V_max[Branchdata[i,2]]+
      lL[i,5]*V_min[Branchdata[i,2]]+lL[i,6]*V_min[Branchdata[i,2]]+lL[i,7]*V_max[Branchdata[i,2]]+lL[i,8]*V_max[Branchdata[i,2]])==0)
      @constraint(model,[i=1:Num_Eline],SR[i]-(lL[i,1]*SR_min[i]+lL[i,2]*SR_max[i]+lL[i,3]*SR_min[i]+lL[i,4]*SR_max[i]+
      lL[i,5]*SR_min[i]+lL[i,6]*SR_max[i]+lL[i,7]*SR_min[i]+lL[i,8]*SR_max[i])==0)
         
   else
         
      Num_pieces=K
      
      #Piecewise quadratic constraints as lambda------------------------------------------------------------------
      KK=Num_pieces+1;
      Xi_V=zeros(Num_bus,KK);
      for i=1:Num_bus
          for m=1:Num_pieces
              Xi_V[i,m]=V_min[i]+((m-1)/Num_pieces).*(V_max[i]-V_min[i]);
          end
      end
      Xi_V[:,KK]=V_max[:];
      
      @variable(model,0<=lq[1:Num_bus,1:KK])
      @variable(model,zq[1:Num_bus,1:Num_pieces],Bin)
      
      @constraint(model,[i=1:Num_bus],sum(zq[i,j] for j=1:Num_pieces)==1)
      @constraint(model,[i=1:Num_bus],sum(lq[i,j] for j=1:KK)==1)
      
      @constraint(model,[i=1:Num_bus],zq[i,1]-lq[i,1]>=0)
      @constraint(model,[i=1:Num_bus],zq[i,KK-1]-lq[i,KK]>=0)
      @constraint(model,[i=1:Num_bus,k=2:KK-1],zq[i,k]+zq[i,k-1]-lq[i,k]>=0)
      
      @constraint(model,[i=1:Num_bus],Vi[i]-sum(lq[i,k]*Xi_V[i,k] for k=1:KK)==0)
      @constraint(model,[i=1:Num_bus],Ui[i]-Vi[i]^2>=0)
      @constraint(model,[i=1:Num_bus],Ui[i]-sum(lq[i,k]*(Xi_V[i,k])^2 for k=1:KK)<=0)
      
      #Strengthening valid inequalities-------------------------------------
      @constraint(model,[i=1:Num_bus],Vi[i]-sum(zq[i,k-1]*Xi_V[i,k] for k=2:KK)<=0)
      @constraint(model,[i=1:Num_bus],Vi[i]-sum(zq[i,k]*Xi_V[i,k] for k=1:KK-1)>=0)
      
      
      #Trilinear constraints as PWM-------------------------------------------------------------------------------
      #----------{Uij=Vi*Vj}--------------
      Vn_min=zeros(Num_bus,Num_pieces);
      Vn_max=zeros(Num_bus,Num_pieces);
      for i=1:Num_bus
          for n=1:Num_pieces
              Vn_min[i,n]=V_min[i]+((n-1)/Num_pieces).*(V_max[i]-V_min[i]);
              Vn_max[i,n]=V_min[i]+((n)/Num_pieces).*(V_max[i]-V_min[i]);
          end
      end
      
      @variable(model,zv[1:Num_bus,1:Num_pieces],Bin)
      @variable(model,Uij[i in 1:Num_Eline])
      @variable(model,Vn[i in 1:Num_bus,1:Num_pieces])
      
      @constraint(model,[i=1:Num_bus],sum(zv[i,j] for j=1:Num_pieces)==1)
      @constraint(model,[i=1:Num_bus,n=1:Num_pieces],Vn[i,n]-Vn_max[i,n].*zv[i,n]<=0)
      @constraint(model,[i=1:Num_bus,n=1:Num_pieces],Vn[i,n]-Vn_min[i,n].*zv[i,n]>=0)
      #@constraint(model,[i=1:Num_bus],Vi[i]-sum(Vn[i,n] for n=1:Num_pieces)==0)
      
      @constraint(model,[i=1:Num_Eline],Uij[i]-sum(Vn_min[Branchdata[i,1],n].*Vn[Branchdata[i,2],n]+Vn[Branchdata[i,1],n].*Vn_min[Branchdata[i,2],n]-Vn_min[Branchdata[i,1],n].*Vn_min[Branchdata[i,2],n].*zv[Branchdata[i,1],n] for n=1:Num_pieces)>=0)
      @constraint(model,[i=1:Num_Eline],Uij[i]-sum(Vn_max[Branchdata[i,1],n].*Vn[Branchdata[i,2],n]+Vn[Branchdata[i,1],n].*Vn_max[Branchdata[i,2],n]-Vn_max[Branchdata[i,1],n].*Vn_max[Branchdata[i,2],n].*zv[Branchdata[i,1],n] for n=1:Num_pieces)>=0)
      @constraint(model,[i=1:Num_Eline],Uij[i]-sum(Vn_max[Branchdata[i,1],n].*Vn[Branchdata[i,2],n]+Vn[Branchdata[i,1],n].*Vn_min[Branchdata[i,2],n]-Vn_max[Branchdata[i,1],n].*Vn_min[Branchdata[i,2],n].*zv[Branchdata[i,1],n] for n=1:Num_pieces)<=0)
      @constraint(model,[i=1:Num_Eline],Uij[i]-sum(Vn_min[Branchdata[i,1],n].*Vn[Branchdata[i,2],n]+Vn[Branchdata[i,1],n].*Vn_max[Branchdata[i,2],n]-Vn_min[Branchdata[i,1],n].*Vn_max[Branchdata[i,2],n].*zv[Branchdata[i,1],n] for n=1:Num_pieces)<=0)
      
      #----------{RE=Uij*CRij}--------------
      CRn_min=zeros(Num_Eline,Num_pieces);
      CRn_max=zeros(Num_Eline,Num_pieces);
      for i=1:Num_Eline
          for n=1:Num_pieces
              CRn_min[i,n]=CR_min[i]+((n-1)/Num_pieces).*(CR_max[i]-CR_min[i]);
              CRn_max[i,n]=CR_min[i]+((n)/Num_pieces).*(CR_max[i]-CR_min[i]);
          end
      end
      
      @variable(model,zcr[1:Num_Eline,1:Num_pieces],Bin)
      @variable(model,CRn[1:Num_Eline,1:Num_pieces])
      
      @constraint(model,[i=1:Num_Eline],sum(zcr[i,j] for j=1:Num_pieces)==1)
      @constraint(model,[i=1:Num_Eline,n=1:Num_pieces],CRn[i,n]-CRn_max[i,n].*zcr[i,n]<=0)
      @constraint(model,[i=1:Num_Eline,n=1:Num_pieces],CRn[i,n]-CRn_min[i,n].*zcr[i,n]>=0)
      @constraint(model,[i=1:Num_Eline],CR[i]-sum(CRn[i,n] for n=1:Num_pieces)==0)
      
      @constraint(model,[i=1:Num_Eline],RE[i]-(V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR[i]+Uij[i]*CR_min[i]-V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR_min[i])>=0)
      @constraint(model,[i=1:Num_Eline],RE[i]-(V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR[i]+Uij[i]*CR_max[i]-V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR_max[i])>=0)
      @constraint(model,[i=1:Num_Eline],RE[i]-(V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR[i]+Uij[i]*CR_min[i]-V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*CR_min[i])<=0)
      @constraint(model,[i=1:Num_Eline],RE[i]-(V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR[i]+Uij[i]*CR_max[i]-V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*CR_max[i])<=0)
      
      #----------{LE=Uij*SRij}--------------
      SRn_min=zeros(Num_Eline,Num_pieces);
      SRn_max=zeros(Num_Eline,Num_pieces);
      for i=1:Num_Eline
          for n=1:Num_pieces
              SRn_min[i,n]=SR_min[i]+((n-1)/Num_pieces).*(SR_max[i]-SR_min[i]);
              SRn_max[i,n]=SR_min[i]+((n)/Num_pieces).*(SR_max[i]-SR_min[i]);
          end
      end
      
      @variable(model,zsr[1:Num_Eline,1:Num_pieces],Bin)
      @variable(model,SRn[1:Num_Eline,1:Num_pieces])
      
      @constraint(model,[i=1:Num_Eline],sum(zsr[i,j] for j=1:Num_pieces)==1)
      @constraint(model,[i=1:Num_Eline,n=1:Num_pieces],SRn[i,n]-SRn_max[i,n].*zsr[i,n]<=0)
      @constraint(model,[i=1:Num_Eline,n=1:Num_pieces],SRn[i,n]-SRn_min[i,n].*zsr[i,n]>=0)
      @constraint(model,[i=1:Num_Eline],SR[i]-sum(SRn[i,n] for n=1:Num_pieces)==0)
      
      @constraint(model,[i=1:Num_Eline],LE[i]-(V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR[i]+Uij[i]*SR_min[i]-V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR_min[i])>=0)
      @constraint(model,[i=1:Num_Eline],LE[i]-(V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR[i]+Uij[i]*SR_max[i]-V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR_max[i])>=0)
      @constraint(model,[i=1:Num_Eline],LE[i]-(V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR[i]+Uij[i]*SR_min[i]-V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]*SR_min[i])<=0)
      @constraint(model,[i=1:Num_Eline],LE[i]-(V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR[i]+Uij[i]*SR_max[i]-V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]*SR_max[i])<=0)
      
   end
   
end

if RT>1
   
   #Candidate line relaxation tightening constraints-------------------------------------------------------------------------
   @variable(model,Theta_min[i+Num_Eline]<=ThetaCL[i in 1:Num_Cline]<=Theta_max[i+Num_Eline])
   @constraint(model,con_ThetaCL[i=1:Num_Cline],ThetaCL[i]-(theta[convert(UInt8,Linedata_candidate[i,1])]-theta[convert(UInt8,Linedata_candidate[i,2])])==0)
   
   #Cosine and Sine constraints---------------
   @variable(model,-1<=CR_CL[1:Num_Cline]<=1)
   @variable(model,-1<=SR_CL[1:Num_Cline]<=1)
     
   Theta_M_CL = max.(abs.(Theta_min[Num_Eline+1:end]), abs.(Theta_max[Num_Eline+1:end]));
   Theta_u_CL=(Theta_min[Num_Eline+1:end].+Theta_max[Num_Eline+1:end])./2;
      
   @constraint(model,con1_CR_CL[i=1:Num_Cline],CR_CL[i]-1+((1-cos(Theta_M_CL[i]))/Theta_M_CL[i]^2)*ThetaCL[i]^2<=0)
   @constraint(model,con2_CR_CL[i=1:Num_Cline],CR_CL[i]-((cos(pi)-cos(-pi))/(2*pi))*(ThetaCL[i]-pi)-cos(pi)>=0)
   @constraint(model,con1_SR_CL[i=1:Num_Cline],SR_CL[i]-cos(Theta_M_CL[i]/2)*(ThetaCL[i]-Theta_M_CL[i]/2)-sin(Theta_M_CL[i]/2)<=0)
   @constraint(model,con2_SR_CL[i=1:Num_Cline],SR_CL[i]-cos(Theta_M_CL[i]/2)*(ThetaCL[i]+Theta_M_CL[i]/2)+sin(Theta_M_CL[i]/2)>=0)
   
   #Trilinear constraints---------------------
   @variable(model,0<=lR_CL[1:Num_Cline,1:8])
   @variable(model,0<=lL_CL[1:Num_Cline,1:8])
   
   CR_min_CL=-ones(Num_Cline);
   CR_max_CL=ones(Num_Cline);
   SR_min_CL=-ones(Num_Cline);
   SR_max_CL=ones(Num_Cline);
   
   @constraint(model,[i=1:Num_Cline],sum(lR_CL[i,j] for j=1:8)==1)
   @constraint(model,[i=1:Num_Cline],sum(lL_CL[i,j] for j=1:8)==1)
   
   @constraint(model,[i=1:Num_Cline],RC[i]-
   (lR_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i]+
   lR_CL[i,3]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,4]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i]+
   lR_CL[i,5]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,6]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i]+
   lR_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i])-1000*(1-XL[i])<=0)
      
   @constraint(model,[i=1:Num_Cline],RC[i]-
   (lR_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i]+
   lR_CL[i,3]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,4]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i]+
   lR_CL[i,5]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,6]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i]+
   lR_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_min_CL[i]+
   lR_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*CR_max_CL[i])+1000*(1-XL[i])>=0)
     
   @constraint(model,[i=1:Num_Cline],Vi[convert(UInt8,Linedata_candidate[i,1])]-(lR_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,3]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,4]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,5]*V_max[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,6]*V_max[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,1])]+lR_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,1])])==0)
   @constraint(model,[i=1:Num_Cline],Vi[convert(UInt8,Linedata_candidate[i,2])]-(lR_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,3]*V_max[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,4]*V_max[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,5]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,6]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,2])]+lR_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,2])])==0)
   @constraint(model,[i=1:Num_Cline],CR_CL[i]-(lR_CL[i,1]*CR_min_CL[i]+lR_CL[i,2]*CR_max_CL[i]+lR_CL[i,3]*CR_min_CL[i]+lR_CL[i,4]*CR_max_CL[i]+lR_CL[i,5]*CR_min_CL[i]+lR_CL[i,6]*CR_max_CL[i]+lR_CL[i,7]*CR_min_CL[i]+lR_CL[i,8]*CR_max_CL[i])==0)  
   
   @constraint(model,[i=1:Num_Cline],LC[i]-
   (lL_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i]+
   lL_CL[i,3]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,4]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i]+
   lL_CL[i,5]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,6]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i]+
   lL_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i])-1000*(1-XL[i])<=0)
   
   @constraint(model,[i=1:Num_Cline],LC[i]-
   (lL_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i]+
   lL_CL[i,3]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,4]*V_min[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i]+
   lL_CL[i,5]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,6]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_min[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i]+
   lL_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_min_CL[i]+
   lL_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,1])]*V_max[convert(UInt8,Linedata_candidate[i,2])]*SR_max_CL[i])+1000*(1-XL[i])>=0)
   
   @constraint(model,[i=1:Num_Cline],Vi[convert(UInt8,Linedata_candidate[i,1])]-(lL_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,3]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,4]*V_min[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,5]*V_max[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,6]*V_max[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,1])]+lL_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,1])])==0)
   @constraint(model,[i=1:Num_Cline],Vi[convert(UInt8,Linedata_candidate[i,2])]-(lL_CL[i,1]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,2]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,3]*V_max[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,4]*V_max[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,5]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,6]*V_min[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,7]*V_max[convert(UInt8,Linedata_candidate[i,2])]+lL_CL[i,8]*V_max[convert(UInt8,Linedata_candidate[i,2])])==0)
   @constraint(model,[i=1:Num_Cline],SR_CL[i]-(lL_CL[i,1]*SR_min_CL[i]+lL_CL[i,2]*SR_max_CL[i]+lL_CL[i,3]*SR_min_CL[i]+lL_CL[i,4]*SR_max_CL[i]+lL_CL[i,5]*SR_min_CL[i]+lL_CL[i,6]*SR_max_CL[i]+lL_CL[i,7]*SR_min_CL[i]+lL_CL[i,8]*SR_max_CL[i])==0)
   
end

#--------------------------------------------------------------------------------------------------------------------------------
@objective(model,Min,sum(Linedata_candidate[i,6]*XL[i] for i in 1:Num_Cline)+(1/100)*sum(Pg[i] for i in 1:Num_gen)+(1/100)*sum(Qi[i] for i in 1:Num_bus));

return model

end


function get_LB(bounds,RT,K)
    model = relax(bounds,RT,K);
    set_optimizer_attribute(model, "CPX_PARAM_EPGAP", 0.0001)#Relative MIP gap tolerance
    #set_optimizer_attribute(model, "CPX_PARAM_PARALLELMODE", 1)#parallel mode switch
    #set_optimizer_attribute(model, "CPX_PARAM_THREADS", 48)#global thread count
    optimize!(model);
    cost=objective_value(model);
    BestBound=objective_bound(model);
    print("XL= ",round.(value.(model[:XL])),"\n\n")
    return BestBound, cost, round.(value.(model[:XL])), value.(model[:Qi])
end
