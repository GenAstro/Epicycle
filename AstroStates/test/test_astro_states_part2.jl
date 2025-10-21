

#  Test an orbit that is a straigth line (pos and vel aligned)
mu = 398600.4418
cartline = CartesianState([7000.0,7000.0,7000.0,7.0,7.0,7.0])
# This should warn and return NaNs
KeplerianState(cartline,mu)
# This should succesfully round trip
to_posvelCartesianState(SphericalRADECState(cartline)) - cartline.posvel
# This should warn and return NaNs
OutGoingAsymptoteState(cartline, mu)





