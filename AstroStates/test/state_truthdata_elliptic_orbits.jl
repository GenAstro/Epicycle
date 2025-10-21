


# GMAT Truth Data
truth_cart = CartesianState([6759.747343616322723, 1115.043329211011041, 1344.722777534846955,
                             -2.660243619064134, 7.541202154282467, 0.640887592324028])

truth_kep = KeplerianState(8000.0, 0.2, deg2rad(12.85), deg2rad(310.0), deg2rad(120.0), deg2rad(300.0))

truth_sphradec = SphericalRADECState(6981.818181818182893, deg2rad(9.366786036853014), deg2rad(11.10476195323573), 
                                     8.022304092374013, deg2rad(109.4308653184559), deg2rad(4.582140635361989))

truth_sphazifpa = SphericalAZIFPAState(6981.818181818182893, deg2rad(9.366786036853014), deg2rad(11.10476195323573), 
                                     8.022304092374013, deg2rad(83.49318158128649), deg2rad(98.94827556462707))

truth_mee = ModifiedEquinoctialState(7680.000000000004547, 0.068404028665134, 0.187938524157182,
                                     0.072384194340782, -0.086264123652664, deg2rad(10.0))

truth_outasymptote = OutGoingAsymptoteState(6400.000000000003, -49.8250551875, deg2rad(250.633213963147),
                                            deg2rad(-11.10476195323574), deg2rad(96.50681841871349), deg2rad(300.0))

# Test case is elliptic so in and out asymptote state are the same (there is no asymptote)
truth_inasymptote = IncomingAsymptoteState(6400.000000000003, -49.8250551875, deg2rad(250.633213963147),
                                            deg2rad(-11.10476195323574), deg2rad(96.50681841871349), deg2rad(300.0));  

truth_modkep = ModifiedKeplerianState(6400.000000000003, 9600.0, deg2rad(12.85),
                                           deg2rad(310.0), deg2rad(120.0), deg2rad(300.0));  

truth_equinoct = EquinoctialState(8000.0 , 0.1879385241571815, 0.0684040286651338, -0.08626412365266437, 0.07238419434078193, deg2rad(28.36066564454829))

truth_altequinoct = AlternateEquinoctialState(8000.0 , 0.1879385241571815, 0.0684040286651338, -0.08572231482903742, 0.07192956275670796, deg2rad(28.36066564454829))

states = Dict(
    :Cartesian => truth_cart,
    :Keplerian => truth_kep,
    :SphericalRADEC => truth_sphradec,
    :SphericalAZIFPAState => truth_sphazifpa,
    :ModifiedEquinoctial => truth_mee,
    :OutGoingAsymptoteState => truth_outasymptote,
    :IncomingAsymptoteState => truth_inasymptote,
    :ModifiedKeplerianState => truth_modkep,
    :EquinoctialState => truth_equinoct,
    :AlternateEquinoctialState => truth_altequinoct
)