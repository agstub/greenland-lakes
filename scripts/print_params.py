# print the auxiliary model parameters

import sys
sys.path.insert(0, '../source')

from params import beta_d,eta_d,t_r,u_d,v_d,H
from metavars import lake_name

def print_params():
    print('lake name: '+lake_name)
    print('\n')
    print('auxiliary model parameters:')
    print('H = '+str(H/1e3)+' km')
    print('beta = '+'{:.2e}'.format(beta_d)+' Pa s m^-1')
    print('eta = '+'{:.2e}'.format(eta_d)+' Pa s')
    print('u = '+str(u_d)+' m yr^-1')
    print('v = '+str(v_d)+' m yr^-1')
    print('relaxation time = '+str(t_r/3.154e7)+' yr')