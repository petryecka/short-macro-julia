


var y           ${y}$        (long_name='output')
    infl        ${\pi}$      (long_name='inflation')
    R           ${R}$        (long_name='Interest Rate')
    g           ${g}$        (long_name='hours')
    z           ${z}$        (long_name='TFP')
    R_obs       ${R_obs}$    (long_name='interest rate observed')
    infl_obs    ${infl_obs}$ (long_name='inflation observed')
    y_obs       ${y_obs}$    (long_name='output observed')
    ;

varexo eps_z ${\varepsilon_z}$ (long_name='TFP shock')
       eps_g ${\varepsilon_g}$ (long_name='government spending shock')
       eps_r ${\varepsilon_r}$ (long_name='monetary policy shock')
    ;
    
parameters 
    tau    ${\tau}$         (long_name='risk aversion')
    kappa  ${\kappa}$       (long_name='phil curve')
    psi1   ${\psi_1}$       (long_name='monetary policy 1')
    psi2   ${\psi_2}$       (long_name='monetary policy 2')
    rhor   ${\rho_r}$       (long_name='mp shock persistence')
    rhog   ${\rho_g}$       (long_name='government shock persistence')
    rhoz   ${\rho_z}$       (long_name='tfp shock persistence')
    ra     ${r^A}$          (long_name='annual interest rate')
    pia    ${pi^A}$         (long_name='annual inflation rate')
    gammaq ${\gamma^Q}$     (long_name='quarterly growth rate')
    sigr   ${\sigma_R}$     (long_name='sd of monetary shock')
    sigg   ${\sigma_g}$     (long_name='sd of gov shock')
    sigz   ${\sigma_z}$     (long_name='sd of tfp')
    ;


tau = 2.83;
kappa = .78;
psi1 = 1.80;
psi2 = .63;
ra = .42;
pia = 3.3;
gammaq = .52;
rhor = .77;
rhog = .98;
rhoz = .88;

model;

[name='Euler equation']

    y = y(+1) - (1/tau)*(R-infl(+1) - z(1)) + g - g(1);

[name='Phillips']

    infl = (1/(1+ra/400))*infl(+1) + kappa*(y - g);

[name='Monetary Policy Rule'] 

    R = rhor*R(-1) + (1-rhor)*psi1*infl + (1-rhor)*psi2*(y-g) + eps_r;

[name='productivity shock']

    z = rhoz*z(-1) + eps_z;

[name='Government shock']

    g = rhog*g(-1) + eps_g;
[name='Output observed']

    y_obs = gammaq + (y - y(-1) + z);

[name='Inflation observed']

    infl_obs = pia + 4*infl;

[name='Interest observed']

    R_obs = pia + ra + 4*gammaq + 4*R;
end;



steady_state_model;
    y=0;
    infl=0;        
    R=0;           
    g=0;         
    z=0;          
    y_obs = gammaq ;
    infl_obs = pia ;
    R_obs = pia + ra + 4*gammaq; 
end;

steady;
check;

//set shock variances
//shocks;
//var eps_z=0.01^2;
//var eps_g=0.01^2;
//var eps_r=0.01^2;
//end;


prior!(:ra, shape=Gamma,mean=0.5,stdev=0.5)

prior!(:pia, shape=Gamma,mean=7,stdev=2)

prior!(:gammaq, shape=Normal, mean=0.4, stdev=0.2)

prior!(:tau, shape=Normal, mean=2, stdev=0.5)

prior!(:kappa, shape=Uniform, mean=0.5, variance=1/12, domain=[0,1])

prior!(:psi1, shape=Normal, mean=1.5, stdev=0.25)

prior!(:psi2, shape=Normal, mean=0.5, stdev=0.25)

prior!(:rhor,shape=Uniform, mean=0.5, variance=1/12, domain=[0,1])

prior!(:rhog, shape=Uniform, mean=0.5, variance=1/12, domain=[0,1])

prior!(:rhoz, shape=Uniform, mean=0.5, variance=1/12, domain=[0,1])



#they need to be declared as standard errors
//prior!(:sigr, shape=InverseGamma1, mean=0.5013256549262002, variance=0.06867258771281654)

//prior!(:sigg, shape=InverseGamma1, mean = 1.2533141373155003, variance=0.4292036732051032)

//prior!(:sigz, shape=InverseGamma1, mean = 0.6266570686577502, variance=0.1073009183012758)

 


estimated_params_init;
  tau, 2.263895645223946;
  kappa, 0.9;
  psi1, 1.932357745633903;
  psi2, 0.465714640873466;
  ra, 0.333453026636814;
  pia, 3.427852381637088;
  gammaq, 0.624767199308368;
  rhor, 0.764476492352786;
  rhog, 0.990359718181570;
  rhoz, 0.917302546854851;
end;

varobs y_obs infl_obs R_obs;

//estimation(datafile='dsge1_data.csv',mh_jscale=0.2,mh_replic=0,mh_nblocks=2);


rwmh_compute!(datafile="dsge1_data.csv", mcmc_replic=100000, mcmc_jscale=0.001, mcmc_chains=1, transformed_parameters = false);

