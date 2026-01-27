def get_model_eqs(sinput):
    eqs_neuron = '''
        # Hodgkin-Huxley for the neuron model (S1L5)
        a_m=-0.182/second/mvolt*(v-U_m)/expm1((U_m-v)/W_m)   : hertz
        b_m=-0.124/second/mvolt*(U_m-v)/expm1((v-U_m)/W_m)   : hertz
        a_h=-0.015/second/mvolt*(U_h-v)/expm1((v-U_h)/W_h)   : hertz
        b_h=-0.015/second/mvolt*(v-U_h)/expm1((U_h-v)/W_h)   : hertz
        m_inf=a_m/(a_m+b_m)                                : 1
        h_inf=a_h/(a_h+b_h)                                : 1
        n_inf=1/(1+exp(-(v-U_n)/W_n))                      : 1
        tau_m=1e-3/(a_m+b_m)/T_adj                         : second
        tau_h=1e-3/(a_h+b_h)/T_adj                         : second
        tau_n=4*ms/(1+exp(-(v+56.56*mvolt)/(44.14*mvolt)))/T_adj      : second
        I_inj=I_dc+I_ramp*t/T_ramp                      : amp/meter**2
        dm/dt=(m_inf-m)/tau_m                           : 1
        dh/dt=(h_inf-h)/tau_h                           : 1
        dn/dt=(n_inf-n)/tau_n                           : 1
    
        # Compute Intracellular Concentrations
        N_n=clip(N0_i+n_Na_n/Lambda_n,0*mmolar,inf*mmolar)  : mmolar
        K_n=clip(K0_i+n_K_n/Lambda_n,0*mmolar,inf*mmolar)   : mmolar
        C_n=clip(C0_i+n_Cl_n/Lambda_n,0*mmolar,inf*mmolar)  : mmolar
    
        # Resolve Nernst Potentials
        V_T=ThermalPotential(T_exp)*volt                : volt
        E_Na=NernstPotential(N_e,N_n,1,T_exp)*volt      : volt
        E_K=NernstPotential(K_e,K_n,1,T_exp)*volt       : volt
        E_Cl=NernstPotential(C_e,C_n,-1,T_exp)*volt     : volt
    
        # Compute Individual fluxes
        I_Na=g_Na*m**3*h*(v-E_Na)                       : amp/meter**2
        I_K=g_K*n*(v-E_K)                               : amp/meter**2
    
        # Leakage components
        I_L_Na=g_L_Na*(v-E_Na)                          : amp/meter**2
        I_L_K=g_L_K*(v-E_K)                             : amp/meter**2
        I_L_Cl=g_L_Cl*(v-E_Cl)                          : amp/meter**2
        I_L=I_L_Na+I_L_K+I_L_Cl                         : amp/meter**2
    
        # Transport mechanisms
        I_NKP=I_NKA*Hill(N_n,zeta_Na,1.5)*Hill(K_e,zeta_K,1)/(1+0.1245*exp(-0.1*v/V_T)-0.0052*exp(-v/V_T)*(1-exp(N_e/67.3/mmolar))) : amp/meter**2
        I_KCC=g_KCC*(E_K-E_Cl)                          : amp/meter**2
        '''

    if sinput=='glu':
        eqs_neuron += '''
                # Synaptic currents
                E_AMPA=V_T*log((N_e+P_K_Na*K_e)/(N_n+P_K_Na*K_n)) : volt
                I_AMPA=G_AMPA*(v-E_AMPA)  : amp/meter**2
                G_AMPA : siemens/meter**2
    
                # ODEs
                dv/dt=(I_inj-I_AMPA-I_Na-I_K-I_L-I_NKP)/c_m       : volt
                dn_Na_n/dt=-S_n*(-pi_AMPA_Na*I_AMPA+I_Na+I_L_Na+3*I_NKP)/F              : mole
                dn_K_n/dt=-S_n*(-pi_AMPA_K*I_AMPA+I_K+I_L_K-2*I_NKP+I_KCC)/F            : mole
                dn_Cl_n/dt=S_n*(I_L_Cl+I_KCC)/F                       : mole
                '''
    else:
        eqs_neuron += '''
                # Synaptic currents
                E_GABA = -V_T*log((C_e+P_HBC_Cl*HBC0_e)/(C_n+P_HBC_Cl*HBC0_i)): volt
                I_GABA = G_GABA*(v-E_GABA): amp/meter**2
                G_GABA: siemens/meter**2

                # ODEs
                dv/dt = (I_inj-pi_GABA_Cl*I_GABA-I_Na-I_K-I_L-I_NKP)/c_m: volt
                dn_Na_n/dt = -S_n*(I_Na+I_L_Na+3*I_NKP)/F: mole
                dn_K_n/dt = -S_n*(I_K+I_L_K-2*I_NKP+I_KCC)/F: mole
                dn_Cl_n/dt = S_n*(pi_GABA_Cl*I_GABA+I_L_Cl+I_KCC)/F: mole
                '''

    eqs_astro = '''    
        # Concentrations
        N_a=clip(N0_a+n_Na_a/Lambda_a,0*mmolar,inf*mmolar)  : mmolar
        K_a=clip(K0_a+n_K_a/Lambda_a,0*mmolar,inf*mmolar)   : mmolar
        C_a=clip(C0_a+n_Cl_a/Lambda_a,0*mmolar,inf*mmolar)  : mmolar
    
        # Define relevant quantities
        E_Na_a=NernstPotential(N_e,N_a,1,T_exp)*volt      : volt 
        E_K_a=NernstPotential(K_e,K_a,1,T_exp)*volt       : volt 
        E_Cl_a=NernstPotential(C_e,C_a,-1,T_exp)*volt     : volt
        E_Glu_a=NernstPotential(G_e,G0_a,-1,T_exp)*volt   : volt
        E_H_a=NernstPotential(H0_e,H0_a,1,T_exp)*volt     : volt
    
        # Kir + XC currents
        I_Kir=g_Kir*(v_a-E_K_a)/(2+exp(1.62*(v_a-E_K_a)/V_T))*Hill(K_e,zeta_Kir,1)   : amp/meter**2
        I_NKP_a=I_NKA_a*Hill(N_a,zeta_Na,1.5)*Hill(K_e,zeta_K,1)/(1+0.1245*exp(-0.1*v_a/V_T)-0.0052*exp(-v_a/V_T)*(1-exp(N_e/67.3/mmolar))) : amp/meter**2
        I_NKCC=g_NKCC*(E_Na_a+E_K_a-2*E_Cl_a)             : amp/meter**2
    
        # Transporter currents
        E_EAAT=(3*E_Na_a+E_H_a-E_K_a-E_Glu_a)/2                 : volt
        E_GAT=(3*E_Na_a+E_Cl_a-V_T*log((GABA0_a/GABA_e)))/2     : volt
        # I_EAAT=g_EAAT*IversonBrackets((G_e-G0_e)/pmolar,1e-9)*(v_a-E_EAAT) : amp/meter**2
        I_EAAT=g_EAAT*int(G_e/G0_e>1.0)*(v_a-E_EAAT) : amp/meter**2
        # I_GAT=g_GAT*IversonBrackets((GABA_e-GABA0_e)/pmolar,1e-9)*(v_a-E_GAT)                           : amp/meter**2
        I_GAT=g_GAT*int(GABA_e/GABA0_e>1.0)*(v_a-E_GAT)                           : amp/meter**2
    
        # Leak current (ClCs + EAATs)
        I_Cl_a=(g_L_Cl_a + g_T_Cl_a*IversonBrackets(abs(I_EAAT/(amp/meter**2)),-1))*(v_a-E_Cl_a)    : amp/meter**2
    
        # GABA-mediated currents
        E_GABA_a=-V_T*log((C_e+P_HBC_Cl*HBC0_e)/(C_a+P_HBC_Cl*HBC0_a)) : volt
        I_GABA_a=g_GABA*r_GABA*(v_a-E_GABA_a)                 : amp/meter**2                
    
        # r.h.s. (astrocyte main (proximal) compartment)    
        dv_a/dt=(-pi_GABA_Cl*I_GABA_a-I_Cl_a-I_Kir-I_NKP_a-I_EAAT-2*I_GAT)/c_m_a : volt
        dn_Na_a/dt=-S_a*(3*I_NKP_a+3*I_EAAT-3*I_GAT-I_NKCC)/F + J_diff_Na*Lambda_a         : mole
        dn_K_a/dt=-S_a*(-I_Kir-2*I_NKP_a-I_EAAT-I_NKCC)/F + J_diff_K*Lambda_a              : mole
        dn_Cl_a/dt=S_a*(-I_GAT+2*I_NKCC+I_Cl_a+pi_GABA_Cl*I_GABA_a)/F + J_diff_Cl*Lambda_a : mole
    
        # Astrocyte-wide GABARs (assuming r_GABA=0 for ambient GABA0_e)              
        dr_GABA/dt = -r_GABA_clip/tau_GABA+J_GABA*(1-r_GABA_clip)*(GABA_e-GABA0_e) : 1        
        r_GABA_clip = clip(r_GABA,0,1)        : 1
        
        # Lumped Buffering (astrocyte)
        J_diff_Na = -D_Na_a*(N_a-N_a_d_clipped) : mmolar/second
        J_diff_K = -D_K_a*(K_a-K_a_d_clipped)   : mmolar/second
        J_diff_Cl = -D_Cl_a*(C_a-C_a_d_clipped) : mmolar/second
        
        # Distal Compartment 
        dN_a_d/dt = -J_diff_Na - J_m_Na : mmolar
        dK_a_d/dt = -J_diff_K  - J_m_K  : mmolar
        dC_a_d/dt = -J_diff_Cl - J_m_Cl : mmolar
        N_a_d_clipped = clip(N_a_d,0*mmolar,inf*mmolar)  : mmolar
        K_a_d_clipped = clip(K_a_d,0*mmolar,inf*mmolar)  : mmolar
        C_a_d_clipped = clip(C_a_d,0*mmolar,inf*mmolar)  : mmolar
        '''

    eqs_ecs = '''
            # ECS-related variations of moles (by diffusion from/to distal compartments)        
            N_e = clip(N0_e+(n_Na_e-n_Na_n-n_Na_a)/Lambda_e, 0*mmolar, inf*mmolar)  : mmolar
            K_e = clip(K0_e+(n_K_e-n_K_n-n_K_a)/Lambda_e, 0*mmolar, inf*mmolar)     : mmolar
            C_e = clip(C0_e+(n_Cl_e-n_Cl_n-n_Cl_a)/Lambda_e, 0*mmolar, inf*mmolar)  : mmolar

            # Lumped buffering (ECS)
            J_diff_Na_e = -D_Na_e*(N_e-N_e_d_clipped)  : mmolar/second
            J_diff_K_e = -D_K_e*(K_e-K_e_d_clipped)    : mmolar/second
            J_diff_Cl_e = -D_Cl_e*(C_e-C_e_d_clipped)  : mmolar/second

            # r.h.s.
            dn_Na_e/dt = J_diff_Na_e*Lambda_e        : mole
            dn_K_e/dt = J_diff_K_e*Lambda_e          : mole
            dn_Cl_e/dt = J_diff_Cl_e*Lambda_e        : mole

            # # Transmembrane fluxes at distal compartments
            J_m_Na = -D_Na_m*(N_e_d_clipped-N_a_d_clipped) : mmolar/second
            J_m_K = -D_K_m*(K_e_d_clipped-K_a_d_clipped)   : mmolar/second
            J_m_Cl = -D_Cl_m*(C_e_d_clipped-C_a_d_clipped) : mmolar/second

            # Distal Compartment (ECS)
            dN_e_d/dt = -J_diff_Na_e + J_m_Na : mmolar
            dK_e_d/dt = -J_diff_K_e + J_m_K   : mmolar
            dC_e_d/dt = -J_diff_Cl_e + J_m_Cl  : mmolar
            N_e_d_clipped = clip(N_e_d,0*mmolar,inf*mmolar)  : mmolar
            K_e_d_clipped = clip(K_e_d,0*mmolar,inf*mmolar)  : mmolar
            C_e_d_clipped = clip(C_e_d,0*mmolar,inf*mmolar)  : mmolar
            '''

    eqs_syn = '''
            # Relevant volumes
            Lambda_e : meter**3 (constant)
            Lambda_s : meter**3 (constant)
            
            # Resting Concentration
            G0_e    : mmolar (constant)
            GABA0_e : mmolar (constant)

            # ECS synaptic concentrations
            G_e = clip(G0_e+n_Glu_e/Lambda_e, 0*mmolar, inf*mmolar)        : mmolar
            GABA_e = clip(GABA0_e+n_GABA_e/Lambda_e, 0*mmolar, inf*mmolar) : mmolar

            # ECS-related Nt clearance by diffusion
            J_diff_Glu_e = -D_Glu_e*(G_e-G0_e)            : mmolar/second
            J_diff_GABA_e = -D_GABA_e*(GABA_e-GABA0_e)            : mmolar/second
            '''
    if sinput=='glu':
        eqs_syn += '''
            # Synaptic concentrations
            G_s = clip(G0_e+n_Glu_s/Lambda_s, 0*mmolar, inf*mmolar)        : mmolar

            # ECS-related variations in neurotransmitter concentrations from/to synapses
            J_diff_Glu_s = -D_Glu_e*(G_s-G_e)             : mmolar/second

            # ODEs        
            dn_Glu_s/dt = Lambda_s*J_diff_Glu_s - R_T*J_r*Lambda_s : mole
            dn_Glu_e/dt = S_a/F*I_EAAT - Lambda_s*J_diff_Glu_s + Lambda_e*J_diff_Glu_e: mole
            dn_GABA_e/dt = S_a/F*I_GAT + Lambda_e*J_diff_GABA_e : mole

            # Synaptic compartment
            J_r     : 1/second
            '''
    else:
        eqs_syn += '''
            # Synaptic concentrations
            GABA_s = clip(GABA0_e+n_GABA_s/Lambda_s, 0*mmolar, inf*mmolar) : mmolar

            # ECS-related variations in neurotransmitter concentrations from/to synapses
            J_diff_GABA_s = -D_GABA_e*(GABA_s-G_e)                : mmolar/second
            
            # ODEs
            dn_Glu_e/dt = Lambda_e*J_diff_Glu_e                      : mole        
            dn_GABA_s/dt = Lambda_s*J_diff_GABA_s - R_T*J_r*Lambda_s : mole
            dn_GABA_e/dt = S_a/F*I_GAT - Lambda_s*J_diff_GABA_s + Lambda_e*J_diff_GABA_e : mole
            
            # Synaptic compartment
            J_r     : 1/second
            '''

    eqs = eqs_neuron + eqs_astro + eqs_ecs + eqs_syn
    return eqs