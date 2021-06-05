Slmax=ELcap*Branchdata[:,6]/baseMVA;
function FBBT!(bounds)
    ep=1e-10;
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

for iter = 1:5
    for i in 1:Num_Eline
    	# PLe_R[i] = GE_tt[i]*Ui[Branchdata[i,2]]+GE_tf[i]*RE[i]-BE_tf[i]*LE[i]
	# PLe_S[i] = GE_ff[i]*Ui[Branchdata[i,1]]+GE_ft[i]*RE[i]+BE_ft[i]*LE[i]
	# QLe_S[i] = -BE_ff[i]*Ui[Branchdata[i,1]]-BE_ft[i]*RE[i]+GE_ft[i]*LE[i]
	# QLe_R[i] = -BE_tt[i]*Ui[Branchdata[i,2]]-BE_tf[i]*RE[i]-GE_tf[i]*LE[i]
	# GE_tf  GE_ft BE_ff BE_tt negative 


        LE_min[i]=max(LE_min[i],((GE_tt[i]*V_min[Branchdata[i,2]]^2+GE_tf[i]*RE_max[i]-Slmax[i])/BE_tf[i])-ep,
				((-GE_ff[i]*V_max[Branchdata[i,1]]^2-GE_ft[i]*RE_min[i]-Slmax[i])/BE_ft[i])-ep,
				(-BE_ff[i]*V_min[Branchdata[i,1]]^2-BE_ft[i]*RE_max[i]-Slmax[i])/(-GE_ft[i])-ep,
				( BE_tt[i]*V_max[Branchdata[i,2]]^2+BE_tf[i]*RE_min[i]-Slmax[i])/(-GE_tf[i])-ep)
        LE_max[i]=min(LE_max[i],((GE_tt[i]*V_max[Branchdata[i,2]]^2+GE_tf[i]*RE_min[i]+Slmax[i])/BE_tf[i])+ep,
				((-GE_ff[i]*V_min[Branchdata[i,1]]^2-GE_ft[i]*RE_max[i]+Slmax[i])/BE_ft[i])+ep,
   				(-BE_ff[i]*V_max[Branchdata[i,1]]^2-BE_ft[i]*RE_min[i]+Slmax[i])/(-GE_ft[i])+ep,
                                ( BE_tt[i]*V_min[Branchdata[i,2]]^2+BE_tf[i]*RE_max[i]+Slmax[i])/(-GE_tf[i])+ep)

        RE_min[i]=max(RE_min[i],((GE_tt[i]*V_min[Branchdata[i,2]]^2-BE_tf[i]*LE_max[i]-Slmax[i])/(-GE_tf[i]))-ep,
				((GE_ff[i]*V_min[Branchdata[i,1]]^2+BE_ft[i]*LE_min[i]-Slmax[i])/(-GE_ft[i]))-ep,
				((-BE_ff[i]*V_min[Branchdata[i,1]]^2+GE_ft[i]*LE_max[i]-Slmax[i])/BE_ft[i])-ep,
				((-BE_tt[i]*V_min[Branchdata[i,2]]^2-GE_tf[i]*LE_min[i]-Slmax[i])/BE_tf[i])-ep)
        RE_max[i]=min(RE_max[i],((GE_tt[i]*V_max[Branchdata[i,2]]^2-BE_tf[i]*LE_min[i]+Slmax[i])/(-GE_tf[i]))+ep,
				((GE_ff[i]*V_max[Branchdata[i,1]]^2+BE_ft[i]*LE_max[i]+Slmax[i])/(-GE_ft[i]))+ep,
				((-BE_ff[i]*V_max[Branchdata[i,1]]^2+GE_ft[i]*LE_min[i]+Slmax[i])/BE_ft[i])+ep,
				((-BE_tt[i]*V_max[Branchdata[i,2]]^2-GE_tf[i]*LE_max[i]+Slmax[i])/BE_tf[i])+ep)


        # PLe_R[i] = GE_tt[i]*Ui[Branchdata[i,2]]+GE_tf[i]*RE[i]-BE_tf[i]*LE[i]
        # PLe_S[i] = GE_ff[i]*Ui[Branchdata[i,1]]+GE_ft[i]*RE[i]+BE_ft[i]*LE[i]
        # QLe_S[i] = -BE_ff[i]*Ui[Branchdata[i,1]]-BE_ft[i]*RE[i]+GE_ft[i]*LE[i]
        # QLe_R[i] = -BE_tt[i]*Ui[Branchdata[i,2]]-BE_tf[i]*RE[i]-GE_tf[i]*LE[i]
        # GE_tf  GE_ft BE_ff BE_tt negative

	Ui1_min = max(0, ((-GE_ft[i]*RE_min[i]-BE_ft[i]*LE_max[i]-Slmax[i])/GE_ff[i]),
		     (BE_ft[i]*RE_min[i]-GE_ft[i]*LE_min[i]-Slmax[i])/(-BE_ff[i]))
        Ui1_max = min(((-GE_ft[i]*RE_max[i]-BE_ft[i]*LE_min[i]+Slmax[i])/GE_ff[i]),
                     (BE_ft[i]*RE_max[i]-GE_ft[i]*LE_max[i]+Slmax[i])/(-BE_ff[i]))

        Ui2_min = max(0, ((-GE_tf[i]*RE_min[i]+BE_tf[i]*LE_min[i]-Slmax[i])/GE_tt[i]),
                     (BE_tf[i]*RE_min[i]+GE_tf[i]*LE_max[i]-Slmax[i])/(-BE_tt[i]))
        Ui2_max = min(((-GE_tf[i]*RE_max[i]+BE_tf[i]*LE_max[i]+Slmax[i])/GE_tt[i]),
                     (BE_tf[i]*RE_max[i]+GE_tf[i]*LE_min[i]+Slmax[i])/(-BE_tt[i]))
           		     
	if Ui1_max <0 || Ui1_max<Ui1_min
	    println("i ",i)
	    println("Ui1_max:", Ui1_max)
	    break
	end
        V_min[Branchdata[i,1]] = max(V_min[Branchdata[i,1]], sqrt(Ui1_min)-ep)
        V_max[Branchdata[i,1]] = min(V_max[Branchdata[i,1]], sqrt(Ui1_max)+ep)
        V_min[Branchdata[i,2]] = max(V_min[Branchdata[i,2]], sqrt(Ui1_min)-ep)
        V_max[Branchdata[i,2]] = min(V_max[Branchdata[i,2]], sqrt(Ui1_max)+ep)


        #CR = cos(Theta)
    	CR_min[i] = max(CR_min[i], min(cos(Theta_min[i]), cos(Theta_max[i]))-ep)
    	CR_max[i] = min(CR_max[i], 1)
    	if Theta_min[i] >=0 ||  Theta_max[i]<=0
            CR_max[i] = min(CR_max[i], max(cos(Theta_min[i]), cos(Theta_max[i]))+ep)
    	end
    	#Theta = acos(CR)
    	if Theta_max[i]<=0
    	    Theta_min[i] = max(Theta_min[i], -acos(CR_min[i])-ep)
            Theta_max[i] = min(Theta_max[i], -acos(CR_max[i])+ep)
    	elseif Theta_min[i] >=0
            Theta_min[i] = max(Theta_min[i],  acos(CR_max[i])-ep)
	    Theta_max[i] = min(Theta_max[i],  acos(CR_min[i])+ep)
    	else
	    Theta_min[i] = max(Theta_min[i], -acos(CR_min[i])-ep)
            Theta_max[i] = min(Theta_max[i],  acos(CR_min[i])+ep)
    	end


        #SR = LE/Vi/Vj
        SR_min[i]=max(SR_min[i], min(LE_min[i]/(V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]), LE_min[i]/(V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]))-ep)
        SR_max[i]=min(SR_max[i], max(LE_max[i]/(V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]), LE_max[i]/(V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]))+ep)

    	#SR = sin(Theta)
    	SR_min[i] = max(SR_min[i], sin(Theta_min[i])-ep)
    	SR_max[i] = min(SR_max[i], sin(Theta_max[i])+ep)
    	#Theta = asin(SR)
    	Theta_min[i] = max(Theta_min[i], asin(SR_min[i])-ep)
    	Theta_max[i] = min(Theta_max[i], asin(SR_max[i])+ep)
    

	#RE = CR*Vi*Vj
    	RE_min[i] = max(RE_min[i], CR_min[i]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]-ep)
    	RE_max[i] = min(RE_max[i], CR_max[i]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]+ep)
    	#LE = SR*Vi*Vj
    	if SR_max[i]<=0
            LE_min[i] = max(LE_min[i], SR_min[i]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]-ep)
            LE_max[i] = min(LE_max[i], SR_max[i]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]+ep)
        elseif SR_min[i]>=0
    	    LE_min[i] = max(LE_min[i], SR_min[i]*V_min[Branchdata[i,1]]*V_min[Branchdata[i,2]]-ep)
    	    LE_max[i] = min(LE_max[i], SR_max[i]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]+ep)
    	else
	    LE_min[i] = max(LE_min[i], SR_min[i]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]-ep)
            LE_max[i] = min(LE_max[i], SR_max[i]*V_max[Branchdata[i,1]]*V_max[Branchdata[i,2]]+ep)
        end

    	#RE^2+LE^2 == UiUj

	if LE_max[i] >= 0 && LE_min[i] <= 0
    	    ss_min = min(RE_min[i]^2, RE_max[i]^2) + 0
	else
	    ss_min = min(RE_min[i]^2, RE_max[i]^2) + min(LE_min[i]^2, LE_max[i]^2)
	end   
    	ss_max = max(RE_min[i]^2, RE_max[i]^2) + max(LE_min[i]^2, LE_max[i]^2)

	Ui1_min = ss_min/V_max[Branchdata[i,2]]/V_max[Branchdata[i,2]]
	Ui1_max = ss_max/V_min[Branchdata[i,2]]/V_min[Branchdata[i,2]]
	Ui2_min = ss_min/V_max[Branchdata[i,1]]/V_max[Branchdata[i,1]]
	Ui2_max = ss_max/V_min[Branchdata[i,1]]/V_min[Branchdata[i,1]]

        V_min[Branchdata[i,1]] = max(V_min[Branchdata[i,1]], sqrt(Ui1_min)-ep)
        V_max[Branchdata[i,1]] = min(V_max[Branchdata[i,1]], sqrt(Ui1_max)+ep)
        V_min[Branchdata[i,2]] = max(V_min[Branchdata[i,2]], sqrt(Ui1_min)-ep)
        V_max[Branchdata[i,2]] = min(V_max[Branchdata[i,2]], sqrt(Ui1_max)+ep)


    	#Theta = atan(LE/RE)
	if RE_min[i] >= 1e-4
    	    at_min = min(LE_min[i]/RE_max[i], LE_min[i]/RE_min[i], LE_max[i]/RE_max[i], LE_max[i]/RE_min[i])
    	    at_max = max(LE_min[i]/RE_max[i], LE_min[i]/RE_min[i], LE_max[i]/RE_max[i], LE_max[i]/RE_min[i])

    	    Theta_min[i] = max(Theta_min[i], atan(at_min)-ep)
    	    Theta_max[i] = min(Theta_max[i], atan(at_max)+ep)
	end
    end
end    
    return bounds;
end


function OBBT!(bounds, UB)
    for iter=1:2
    	FBBT!(bounds);

	println("FBBT bounds.V_min  bounds.V_max: ", bounds.V_min, bounds.V_max)
	println("FBBT bounds.Theta_min  bounds.Theta_max: ", bounds.Theta_min, bounds.Theta_max)
	println("FBBT bounds.RE_min  bounds.RE_max: ", bounds.RE_min, bounds.RE_max)
	println("FBBT bounds.LE_min  bounds.LE_max: ", bounds.LE_min, bounds.LE_max)
       	println("FBBT bounds.CR_min  bounds.CR_max: ", bounds.CR_min, bounds.CR_max)
        println("FBBT bounds.SR_min  bounds.SR_max: ", bounds.SR_min, bounds.SR_max)

	@everywhere epsilon = 1e-4
        @eval @everywhere bounds = $bounds
        @eval @everywhere UB = $UB
        @everywhere function local_V(i)
                ob = relax(bounds,RT,K);
                @constraint(ob,sum(Linedata_candidate[j,6]*ob[:XL][j] for j in 1:Num_Cline)+(1/100)*sum(ob[:Pg][j] for j in 1:Num_gen)<=UB);
                set_optimizer_attribute(ob, "CPX_PARAM_TILIM", 120)#Time limit
                set_optimizer_attribute(ob, "CPX_PARAM_SCRIND", false)#NO solver print
                set_optimizer_attribute(ob, "CPX_PARAM_PARALLELMODE", 0)#parallel mode switch

                @objective(ob, Min, ob[:Vi][i])
                optimize!(ob);
                min_value = max(objective_bound(ob) - epsilon, bounds.V_min[i])

                println(i, "  max")
                @objective(ob, Max, ob[:Vi][i])
                optimize!(ob);
                max_value = min(objective_bound(ob) + epsilon, bounds.V_max[i])

                @objective(ob, Min, ob[:Ui][i])
                optimize!(ob);
                min_value = max(sqrt(objective_bound(ob)) - epsilon, min_value)

                @objective(ob, Max, ob[:Ui][i])
                optimize!(ob);
                max_value = min(sqrt(objective_bound(ob)) + epsilon, max_value)
		return [min_value, max_value]
        end
        res = pmap(local_V, 1:Num_bus)
        red = ones(Num_bus)
        for i=1:Num_bus
            if (bounds.V_max[i]-bounds.V_min[i])>=1e-4
                min_value  =  res[i][1]
                max_value  =  res[i][2]
                red[i] = (max_value-min_value)/(bounds.V_max[i]-bounds.V_min[i])
                bounds.V_min[i] = min_value
                bounds.V_max[i] = max_value
            end
        end
        println("OBBT bounds.V_min  bounds.V_max: ", mean(red))
	println(bounds.V_min, bounds.V_max)


        @eval @everywhere bounds = $bounds
        @everywhere function local_Theta(i)
                ob = relax(bounds,RT,K);
                @constraint(ob,sum(Linedata_candidate[j,6]*ob[:XL][j] for j in 1:Num_Cline)+(1/100)*sum(ob[:Pg][j] for j in 1:Num_gen)<=UB);
                set_optimizer_attribute(ob, "CPX_PARAM_TILIM", 120)#Time limit
                set_optimizer_attribute(ob, "CPX_PARAM_SCRIND", false)#NO solver print
                set_optimizer_attribute(ob, "CPX_PARAM_PARALLELMODE", 0)#parallel mode switch

                @objective(ob, Min, ob[:Theta][i])
                optimize!(ob);
                min_value = max(objective_bound(ob) - epsilon, bounds.Theta_min[i])

                @objective(ob, Max, ob[:Theta][i])
                optimize!(ob);
                max_value = min(objective_bound(ob) + epsilon, bounds.Theta_max[i])
		return [min_value, max_value]
        end
        res = pmap(local_Theta, 1:Num_Eline)
        red = ones(Num_Eline)
        for i=1:Num_Eline
	    if (bounds.Theta_max[i]- bounds.Theta_min[i])>=1e-4
                min_value  =  res[i][1]
                max_value  =  res[i][2]
                red[i] = (max_value-min_value)/(bounds.Theta_max[i]- bounds.Theta_min[i])
                bounds.Theta_min[i] = min_value
                bounds.Theta_max[i] = max_value
            end
        end
	println("OBBT bounds.Theta_min  bounds.Theta_max: ",  mean(red))
	println(bounds.Theta_min, bounds.Theta_max)


        @eval @everywhere bounds = $bounds
        @everywhere function local_RE(i)
                ob = relax(bounds,RT,K);
                @constraint(ob,sum(Linedata_candidate[j,6]*ob[:XL][j] for j in 1:Num_Cline)+(1/100)*sum(ob[:Pg][j] for j in 1:Num_gen)<=UB);
                set_optimizer_attribute(ob, "CPX_PARAM_TILIM", 120)#Time limit
                set_optimizer_attribute(ob, "CPX_PARAM_SCRIND", false)#NO solver print
                set_optimizer_attribute(ob, "CPX_PARAM_PARALLELMODE", 0)#parallel mode switch

                @objective(ob, Min, ob[:RE][i])
                optimize!(ob);
                min_value = max(objective_bound(ob) - epsilon, bounds.RE_min[i])

                @objective(ob, Max, ob[:RE][i])
                optimize!(ob);
                max_value = min(objective_bound(ob) + epsilon, bounds.RE_max[i])
                return [min_value, max_value]
        end
        res = pmap(local_RE, 1:Num_Eline)
        red = ones(Num_Eline)
        for i=1:Num_Eline
            if (bounds.RE_max[i]- bounds.RE_min[i])>=1e-4
                min_value  =  res[i][1]
                max_value  =  res[i][2]
                red[i] = (max_value-min_value)/(bounds.RE_max[i]- bounds.RE_min[i])
                bounds.RE_min[i] = min_value
                bounds.RE_max[i] = max_value
            end
        end
        println("OBBT bounds.RE_min  bounds.RE_max: ",  mean(red))
        println(bounds.RE_min, bounds.RE_max)


        @eval @everywhere bounds = $bounds
        @everywhere function local_LE(i)
                ob = relax(bounds,RT,K);
                @constraint(ob,sum(Linedata_candidate[j,6]*ob[:XL][j] for j in 1:Num_Cline)+(1/100)*sum(ob[:Pg][j] for j in 1:Num_gen)<=UB);
                set_optimizer_attribute(ob, "CPX_PARAM_TILIM", 120)#Time limit
                set_optimizer_attribute(ob, "CPX_PARAM_SCRIND", false)#NO solver print
                set_optimizer_attribute(ob, "CPX_PARAM_PARALLELMODE", 0)#parallel mode switch

                @objective(ob, Min, ob[:LE][i])
                optimize!(ob);
                min_value = max(objective_bound(ob) - epsilon, bounds.LE_min[i])

                @objective(ob, Max, ob[:LE][i])
                optimize!(ob);
                max_value = min(objective_bound(ob) + epsilon, bounds.LE_max[i])
                return [min_value, max_value]
        end
        res = pmap(local_LE, 1:Num_Eline)
        red = ones(Num_Eline)
        for i=1:Num_Eline
            if (bounds.LE_max[i]- bounds.LE_min[i])>=1e-4
                min_value  =  res[i][1]
                max_value  =  res[i][2]
                red[i] = (max_value-min_value)/(bounds.LE_max[i]- bounds.LE_min[i])
                bounds.LE_min[i] = min_value
                bounds.LE_max[i] = max_value
            end
        end
        println("OBBT bounds.LE_min  bounds.LE_max: ",  mean(red))
        println(bounds.LE_min, bounds.LE_max)



        @eval @everywhere bounds = $bounds
        @everywhere function local_CR(i)
                ob = relax(bounds,RT,K);
                @constraint(ob,sum(Linedata_candidate[j,6]*ob[:XL][j] for j in 1:Num_Cline)+(1/100)*sum(ob[:Pg][j] for j in 1:Num_gen)<=UB);
                set_optimizer_attribute(ob, "CPX_PARAM_TILIM", 120)#Time limit
                set_optimizer_attribute(ob, "CPX_PARAM_SCRIND", false)#NO solver print
                set_optimizer_attribute(ob, "CPX_PARAM_PARALLELMODE", 0)#parallel mode switch

                @objective(ob, Min, ob[:CR][i])
                optimize!(ob);
                min_value = max(objective_bound(ob) - epsilon, bounds.CR_min[i])

                @objective(ob, Max, ob[:CR][i])
                optimize!(ob);
                max_value = min(objective_bound(ob) + epsilon, bounds.CR_max[i])
                return [min_value, max_value]
        end
        res = pmap(local_CR, 1:Num_Eline)
        red = ones(Num_Eline)
        for i=1:Num_Eline
            if (bounds.CR_max[i]- bounds.CR_min[i])>=1e-4
                min_value  =  res[i][1]
                max_value  =  res[i][2]
                red[i] = (max_value-min_value)/(bounds.CR_max[i]- bounds.CR_min[i])
                bounds.CR_min[i] = min_value
                bounds.CR_max[i] = max_value
            end
        end
        println("OBBT bounds.CR_min  bounds.CR_max: ",  mean(red))
        println(bounds.CR_min, bounds.CR_max)


        @eval @everywhere bounds = $bounds
        @everywhere function local_SR(i)
                ob = relax(bounds,RT,K);
                @constraint(ob,sum(Linedata_candidate[j,6]*ob[:XL][j] for j in 1:Num_Cline)+(1/100)*sum(ob[:Pg][j] for j in 1:Num_gen)<=UB);
                set_optimizer_attribute(ob, "CPX_PARAM_TILIM", 120)#Time limit
                set_optimizer_attribute(ob, "CPX_PARAM_SCRIND", false)#NO solver print
                set_optimizer_attribute(ob, "CPX_PARAM_PARALLELMODE", 0)#parallel mode switch

                @objective(ob, Min, ob[:SR][i])
                optimize!(ob);
                min_value = max(objective_bound(ob) - epsilon, bounds.SR_min[i])

                @objective(ob, Max, ob[:SR][i])
                optimize!(ob);
                max_value = min(objective_bound(ob) + epsilon, bounds.SR_max[i])
                return [min_value, max_value]
        end
        res = pmap(local_SR, 1:Num_Eline)
        red = ones(Num_Eline)
        for i=1:Num_Eline
            if (bounds.SR_max[i]- bounds.SR_min[i])>=1e-4
                min_value  =  res[i][1]
                max_value  =  res[i][2]
                red[i] = (max_value-min_value)/(bounds.SR_max[i]- bounds.SR_min[i])
                bounds.SR_min[i] = min_value
                bounds.SR_max[i] = max_value
            end
        end
        println("OBBT bounds.SR_min  bounds.SR_max: ",  mean(red))
        println(bounds.SR_min, bounds.SR_max)

    end
    return bounds;
end






