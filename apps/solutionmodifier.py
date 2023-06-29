import cantera as ct

### Constants
Avogadro = 6.0221409e+23; #[1/mol];
rho_soot = 1800; #[kg/m3];
kB = 1.3806488e-23; #[m2.kg/s2-K]
Pi = 3.14159265359;
C = 2.2 * (8 * Pi * kB) ** 0.5;

def text_num(name):
    head = name.rstrip('0123456789')
    tail = name[len(head):]
    return head, tail


def k_f1_dict(PAH, **kwargs):
    PAH_n_C = PAH.composition["C"];
    A = 9.8e1 * PAH_n_C; #[mol/m3-s] 
    b = 1.84;
    Ea = 62886.6; #[j/mol];
    return {"rate-constant": {"A": A*1e3, "b": b, "Ea": Ea*1e3}};

#----------------------------------------------------------------------------
### Reaction rate constant dict constructor functions
def k_r1_dict(**kwargs):
#     Ref: Semenikhin et al.(2017); Rate constants for H abstraction from benzo(a)pyrene and chrysene: a theoretical study
    A = 1.60e-2;
    b = 2.63;
    Ea = 17837.4; #[j/mol];
    return {"rate-constant": {"A": A*1e3, "b": b, "Ea": Ea*1e3}};

def k_f2_dict(**kwargs):
    #   Ref: Hadrings (2005); Predictive Theory for Hydrogen Atom-Hydrocarbon Radical Association Kinetics
    A = 8.08e-11 * 1e-6 * Avogadro;
    b = 0.13;
    Ea = 0.0;
    return {"rate-constant": {"A": A*1e3, "b": b, "Ea": Ea*1e3}};

def k_f3_dict(PAH_i, PAH_j, **kwargs):
    PAH_mass_i = PAH_i.molecular_weight / (Avogadro*1e3); #kg
    PAH_mass_j = PAH_j.molecular_weight / (Avogadro*1e3); #kg
    d_mi = (6*PAH_mass_i/(Pi*rho_soot))**(1.0/3.0);
    d_mj = (6*PAH_mass_j/(Pi*rho_soot))**(1.0/3.0);
    mu = PAH_mass_i * PAH_mass_j / (PAH_mass_i + PAH_mass_j);
    A = Avogadro * C * (1.0 / mu)**0.5 * (d_mi + d_mj) ** 2.0;
    b = 0.5;
    Ea = 0.0;
    return {"rate-constant": {"A": A*1e3, "b": b, "Ea": Ea*1e3}};


def k_r3_dict(**kwargs):
    A = 1.3e83;
    b = -20.162;
    Ea = 401287.44;
    return {"rate-constant": {"A": A*1e3, "b": b, "Ea": Ea*1e3}};
#----------------------------------------------------------------------------

def add_e_bridge(gas, precursor_list, n_steps = 4):
    e_bridge_species = [];
    e_bridge_reactions = [];

    for precursor_name in precursor_list:
        precursor = gas.species(precursor_name);
        for polymer_index in range(n_steps):
            if polymer_index == 0:
                PAH = precursor;
            else:
                PAH = polymer;
            ## New radical species
            radical_dict = {
                "name":  PAH.name + "e-",
                "composition": {"C": PAH.composition["C"], "H": (PAH.composition["H"]-1)},
                "charge" : PAH.charge,
                "size" : PAH.size,
                "thermo" : None,
                "transport" : None,
            }
            radical = ct.Species(**radical_dict);
            radical.thermo = PAH.thermo;
            radical.transport = PAH.transport;
            e_bridge_species.append(radical);
            ## PAH dehydrogenation forward reaction
            reaction_rate = ct.ReactionRate.from_dict(k_f1_dict(PAH))
            reaction_dict = {
                "reactants" : {PAH.name: 1.0, "H": 1.0},
                "products"  : {radical.name: 1.0, "H2": 1.0},
                "rate"      : reaction_rate,
            }
            reaction = ct.Reaction(**reaction_dict)
            reaction.reversible = False;
            e_bridge_reactions.append(reaction);
            ## PAH dehydrogenation reverse reaction
            reaction_rate = ct.ReactionRate.from_dict(k_r1_dict())
            reaction_dict = {
                "reactants" : {radical.name: 1.0, "H2": 1.0},
                "products"  : {PAH.name: 1.0, "H": 1.0},
                "rate"      : reaction_rate,
            }
            reaction = ct.Reaction(**reaction_dict)
            reaction.reversible = False;
            e_bridge_reactions.append(reaction);
            ## Radical hydrogenation forward reaction
            reaction_rate = ct.ReactionRate.from_dict(k_f2_dict())
            reaction_dict = {
                "reactants" : {radical.name: 1.0, "H": 1.0},
                "products"  : {PAH.name: 1.0},
                "rate"      : reaction_rate,
            }
            reaction = ct.Reaction(**reaction_dict)
            reaction.reversible = False;
            e_bridge_reactions.append(reaction);
            if (polymer_index < (n_steps-1)):
                ## New Polymer
                polymer_dict = {
                    "name":  f"{precursor.name}e{(polymer_index+2)}",
                    "composition": {"C": (radical.composition["C"]+precursor.composition["C"]),
                                    "H": (radical.composition["H"]+precursor.composition["H"]-1)},
                    "charge" : PAH.charge,
                    "size" : PAH.size,
                    "thermo" : None,
                    "transport" : None,
                };
                polymer = ct.Species(**polymer_dict);
                polymer.thermo = PAH.thermo;
                polymer.transport = PAH.transport;
                e_bridge_species.append(polymer);
                ## Polymerization forward reaction
                reaction_rate = ct.ReactionRate.from_dict(k_f3_dict(radical, precursor))
                reaction_dict = {
                    "reactants" : {radical.name: 1.0, precursor.name: 1.0},
                    "products"  : {polymer.name: 1.0, "H": 1.0},
                    "rate"      : reaction_rate,
                }
                reaction = ct.Reaction(**reaction_dict)
                reaction.reversible = False;
                e_bridge_reactions.append(reaction);
                ## Polymerization reverse reaction
                reaction_rate = ct.ReactionRate.from_dict(k_r3_dict())
                reaction_dict = {
                    "reactants" : {polymer.name: 1.0, "H": 1.0},
                    "products"  : {radical.name: 1.0, precursor.name: 1.0},
                    "rate"      : reaction_rate,
                }
                reaction = ct.Reaction(**reaction_dict)
                reaction.reversible = False;
                e_bridge_reactions.append(reaction);
    solution_dict = dict(
        thermo    = 'IdealGas',
        kinetics  = 'GasKinetics',
        transport_model = 'Mix',
        species   =  gas.species() + e_bridge_species,
        reactions =  gas.reactions() + e_bridge_reactions,
    )
    gas_e_bridge= ct.Solution(**solution_dict);
    return gas_e_bridge, e_bridge_species, e_bridge_reactions;