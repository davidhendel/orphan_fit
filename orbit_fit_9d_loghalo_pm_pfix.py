#emcee orbit fit

import sys
import os
#import cld_fnc as cf
from astropy.io import ascii
from astropy.coordinates import ICRS, Galactic, SkyCoord, Distance
from astropy import units as u
import astropy.coordinates as coord
import gala.coordinates as gc
from gala.units import galactic
from scipy.interpolate import interp1d
import scipy
import gala
import astropy.constants as ac
from scipy.interpolate import interp1d
import emcee
import corner
import gala.coordinates as gc
import matplotlib.cm as cm
import numpy as np
from numpy import *

sys.path.append(os.getcwd()) 

pm_prior_on=True


if len(sys.argv) > 2:
	simname = str(sys.argv[1])
	table = np.load('simtable_'+simname+'.npy')
	v_cArray, r_hArray = np.load('vc_rh_arrays_207_log.npy')
	rh_interp = interp1d(v_cArray,r_hArray)
	data_pm_l_cosb = float(sys.argv[2])
	data_pm_b = float(sys.argv[3])
	if pm_prior_on:
		outname = 'samples_log_'+simname+'_pfix_pm.npy'
	else: outname = 'samples_log_'+simname+'_pfix.npy'

else:
	#names should be jill or jill_no19
	name = str(sys.argv[1])
	table = np.load('smhash_rrly_datatable_'+name+'.npy')
	data_pm_l_cosb =  0.211
	data_pm_b = -0.774
	v_cArray, r_hArray = np.load('vc_rh_arrays_log.npy')
	rh_interp = interp1d(v_cArray,r_hArray)
	if pm_prior_on:
		outname = 'samples_log_'+name+'_pfix_pm.npy'
	else: outname = 'samples_log_'+name+'_pfix.npy'


def distance_simple(apmag, apmag_err, period, z, z_slope = 0.184):

	z_err = 0.15
	absmag = -2.276*np.log10(period) + z_slope*z - 0.786
	sigma_zterm = z_slope*z * np.sqrt((0.004/z_slope)**2 + (z_err/z)**2)
	use_zp_err = 1.
	absmag_err = np.sqrt(sigma_zterm**2 + (.007*(use_zp_err))**2 + 0.021**2 + 0.035**2)
	mu = apmag - absmag
	sigma_mu = np.sqrt(absmag_err**2 + apmag_err**2)
	distance = 10**((mu/5.) + 1)
	distance_err = 0.461*distance*sigma_mu

	return distance, distance_err

def distance_modulus(apmag, apmag_err, period, z, z_slope = 0.184):

	z_err = 0.15
	absmag = -2.276*np.log10(period) + z_slope*z - 0.786
	sigma_zterm = z_slope*z * np.sqrt((0.004/z_slope)**2 + (z_err/z)**2)
	use_zp_err = 1.
	absmag_err = np.sqrt(sigma_zterm**2 + (.007*(use_zp_err))**2 + 0.021**2 + 0.035**2)
	mu = apmag - absmag
	sigma_mu = np.sqrt(absmag_err**2 + apmag_err**2)
	distance = 10**((mu/5.) + 1)
	distance_err = 0.461*distance*sigma_mu

	return mu, sigma_mu

#make Newberg potential
pot = gala.potential.CompositePotential()
pot['bulge'] = gala.potential.HernquistPotential(m=3.4e10, c=0.7, units=(u.kpc,u.Myr,u.Msun,u.rad))
pot['disk']  = gala.potential.MiyamotoNagaiPotential(m=1.0e11,a=6.5,b=0.26, units=(u.kpc,u.Myr,u.Msun,u.rad))
pot['halo']  = gala.potential.LogarithmicPotential(v_c=(73.*np.sqrt(2.)*u.km/u.second).to(u.kpc/u.Myr).value, r_h=12.,q1=1., q2=1., q3=1., units=(u.kpc,u.Myr,u.Msun,u.rad))

#if 0:
#	print "dkafsdlj"
	#####################
	# v_c - r_h lookup table

	#def get_rotv(r_h, v_c):
	#	pot['halo']  = gala.potential.LogarithmicPotential(v_c=(v_c*np.sqrt(2.)*u.km/u.second).to(u.kpc/u.Myr).value, r_h=r_h, q1=1., q2=1., q3=1., units=(u.kpc,u.Myr,u.Msun,u.rad))
	#	return (pot.circular_velocity([-8,0,0]).to(u.km/u.s) - 220.*u.km/u.s).value
	#
	#v_cArray = np.arange(70,200,.1)
	#r_hArray = np.zeros(len(v_cArray))
	#for i in np.arange(len(v_cArray)):
	#	try: r_hArray[i] = scipy.optimize.bisect(get_rotv, 0.1, 100., args=(v_cArray[i]))
	#	except: continue
	#
	#np.save('vc_rh_arrays.npy', (v_cArray, r_hArray))
	#
	#plt.plot(v_cArray,r_hArray)
	#plt.xlabel('$v_c$', size='x-large')
	#plt.ylabel('$r_h$', size='x-large')
	#plt.title('$r_h\ s.t. \ v_{rot} |_{R=8kpc} = 220 km/s$', size='x-large')
	#
	#pos = np.zeros((3,1000)) * u.kpc
	#pos[0] = np.linspace(0.1,30.1,pos.shape[1]) * u.kpc
	#norm = mpl.colors.Normalize(vmin=70., vmax=141.)
	#m = cm.ScalarMappable(norm=norm, cmap=cm.rainbow)
	#
	#for i in np.arange(len(v_cArray)):
	#	pot['halo']  = gala.potential.LogarithmicPotential(v_c=(v_cArray[i]*np.sqrt(2.)*u.km/u.second).to(u.kpc/u.Myr).value, r_h=r_hArray[i], q1=1., q2=1., q3=1., units=(u.kpc,u.Myr,u.Msun,u.rad))
	#	vrots = pot.circular_velocity(pos)
	#	plt.plot(pos[0],vrots.to(u.km/u.s), c=m.to_rgba(v_cArray[i]))
	#
	#plt.xlabel('R [kpc]')
	#plt.ylabel('Circular velocity [km/s]')
	#plt.ylim([200,240])
	#plt.xlim([0,30])

	#####################


#Sohn IC
sohn_ic = SkyCoord(ra='10h03m48.9s', dec = '+29d06m12.8s')

#newberg ICs
(l,b,R)=(218., 53.5, 28.6)
newberg_ic = SkyCoord(l = l, b= b, distance=R, unit=(u.degree, u.degree, u.kpc),frame='galactic', galcen_distance=8.0*u.kpc)
rr_c = SkyCoord(ra=table['R.A.'], dec = table['Decl.'], distance=table['Helio. Distance'], unit=(u.degree,u.degree,u.kpc),galcen_distance=8.0*u.kpc)
rr_cg = rr_c.galactic
orp = rr_c.transform_to(gc.Orphan)
orp.Lambda.wrap_angle=180*u.degree
#ttt = table['v_gsr']*(u.km/u.s)
orp_vr = table['v_helio']*(u.km/u.s)#gc.vgsr_to_vhel(orp,ttt)

newberg_init = [newberg_ic.cartesian.x.value-8., 
		newberg_ic.cartesian.y.value, 
		newberg_ic.cartesian.z.value, 
		(-156/977.7922216731282), 
		(79/977.7922216731282), 
		(107/977.7922216731282)]

c = coord.Galactocentric([newberg_ic.cartesian.x - 8*u.kpc]*u.kpc, 
		[newberg_ic.cartesian.y]*u.kpc, 
		[newberg_ic.cartesian.z]*u.kpc, galcen_distance=8.0* u.kpc)

newberg_vxyz = [[-156.], [79.], [107.]]*u.km/u.s
c_gc = c.transform_to(coord.Galactic)
c_icrs = c.transform_to(coord.ICRS)
pm_l, pm_b, vrad = gc.vgal_to_hel(c_gc, newberg_vxyz)
newberg_reverse_orbit = pot.integrate_orbit(array(newberg_init), dt=-1, nsteps=500)
newberg_orbit = pot.integrate_orbit(newberg_reverse_orbit[-1], dt= 1, nsteps=1000)
empty_galframe = coord.Galactocentric(galcen_distance=8.*u.kpc)
newberg_orbitg = newberg_orbit.to_coord_frame(coord.Galactic, galactocentric_frame=empty_galframe)

newberg_interp_distance = interp1d(newberg_orbitg.l, newberg_orbitg.distance)
newberg_interp_x = interp1d(newberg_orbitg.l, newberg_orbit.pos.x)
newberg_interp_y = interp1d(newberg_orbitg.l, newberg_orbit.pos.y)
newberg_interp_z = interp1d(newberg_orbitg.l, newberg_orbit.pos.z)
newberg_interp_vx = interp1d(newberg_orbitg.l, newberg_orbit.vel.d_x.to(u.km/u.s))
newberg_interp_vy = interp1d(newberg_orbitg.l, newberg_orbit.vel.d_y.to(u.km/u.s))
newberg_interp_vz = interp1d(newberg_orbitg.l, newberg_orbit.vel.d_z.to(u.km/u.s))
newberg_interp_pm_l_cosb = interp1d(newberg_orbitg.l, newberg_orbitg.pm_l_cosb)
newberg_interp_pm_b = interp1d(newberg_orbitg.l, newberg_orbitg.pm_b)
newberg_interp_vrad = interp1d(newberg_orbitg.l, newberg_orbitg.radial_velocity)

newberg_orbit_icrs =  newberg_orbitg.transform_to(coord.ICRS)
newberg_interp_pm_ra_cosdec = interp1d(newberg_orbitg.l, newberg_orbit_icrs.pm_ra_cosdec)
newberg_interp_pm_dec = interp1d(newberg_orbitg.l, newberg_orbit_icrs.pm_dec)
newberg_interp_vrad_icrs = interp1d(newberg_orbitg.l, newberg_orbit_icrs.radial_velocity)

sohn_init = [newberg_interp_x(sohn_ic.galactic.l),
		newberg_interp_y(sohn_ic.galactic.l),
		newberg_interp_z(sohn_ic.galactic.l),
		newberg_interp_vx(sohn_ic.galactic.l),
		newberg_interp_vy(sohn_ic.galactic.l),
		newberg_interp_vz(sohn_ic.galactic.l)]

sohn_c = SkyCoord(ra='10h03m48.9s', dec = '+29d06m12.8s', distance = newberg_interp_distance(sohn_ic.galactic.l),unit=(u.degree, u.degree, u.kpc),frame='icrs', galcen_distance=8.0*u.kpc)

sohn_icrs = coord.ICRS(ra=150.94036304*u.deg,  dec=29.10690538*u.deg, pm_ra_cosdec=-.74*u.mas/u.yr, pm_dec=-.31*u.mas/u.yr)

sohn_cg = c.transform_to(coord.Galactic)
sohn_gc = SkyCoord(l=199.77961451,  b=53.18030411,  distance=33.98261686, unit=(u.degree, u.degree, u.kpc),frame='galactic', galcen_distance=8.0*u.kpc)
sohn_vxyz = ([newberg_interp_vx(sohn_ic.galactic.l), newberg_interp_vy(sohn_ic.galactic.l), newberg_interp_vz(sohn_ic.galactic.l)]*u.kpc/u.Myr).to(u.km/u.s)
sohn_vxyz = [[sohn_vxyz[0].value],[sohn_vxyz[1].value],[sohn_vxyz[2].value]]*u.km/u.s
sohn_reverse_orbit = pot.integrate_orbit(array(sohn_init), dt=-1, nsteps=500)
sohn_orbit = pot.integrate_orbit(sohn_reverse_orbit[-1], dt= 1, nsteps=1000)
#sohn_orbitg, sohn_vels = sohn_orbit.to_coord_frame(coord.Galactic)

soh_icrs = coord.ICRS(ra='10h03m48.9s', dec = '+29d06m12.8s', pm_ra_cosdec=-0.74*u.mas/u.yr, pm_dec=-.31*u.mas/u.yr)
sohn_gal = sohn_icrs.transform_to(coord.Galactic)

if 'avg. mag ' in table.dtype.names:
	DMs, DMs_err = distance_modulus(table['avg. mag'], table['mag err'], table['Period'], table['[Fe/H]'])
else:
	DMs = 5.*np.log10(table['Helio. Distance']*1000.)-5.
	DMs_err = np.zeros(len(DMs)) + 0.055

# p = l,b,distance,pm_l,pm_b,v_r

def log_likelihood(p, nsteps=500, l0 = sohn_gal.l.value, ls = orp.Lambda, bs= orp.Beta, DMs = DMs, DMs_err = DMs_err, vr=orp_vr, e_decs = 0., e_vr = table['v_helio err']):

	l = l0
	b, DM, pm_l, pm_b, vrad, sb, sv, sdm, vc = p
	try: rh_interp(vc)
	except ValueError:
		return -inf

	pot['halo']  = gala.potential.LogarithmicPotential(v_c=(vc*np.sqrt(2.)*u.km/u.second).to(u.kpc/u.Myr).value, r_h=rh_interp(vc), q1=1., q2=1., q3=1., units=(u.kpc,u.Myr,u.Msun,u.rad))
	
	distance = 10.**((DM/5.)+1.)/1000.

	#print DM, distance
	try: new_cg = coord.Galactic(l=l*u.deg,b=b*u.deg,distance=distance*u.kpc,pm_l_cosb=pm_l*u.mas/u.yr,pm_b=pm_b*u.mas/u.yr,radial_velocity=vrad*u.km/u.s)
	except:
		return -inf

	new_cg_galactocentric = new_cg.transform_to(coord.Galactocentric(galcen_distance=8.*u.kpc))

	new_init = [new_cg_galactocentric.x.value, 
			new_cg_galactocentric.y.value, 
			new_cg_galactocentric.z.value, 
			(new_cg_galactocentric.v_x.value/977.7922216731282), 
			(new_cg_galactocentric.v_y.value/977.7922216731282), 
			(new_cg_galactocentric.v_z.value/977.7922216731282)]

	reverse_orbit = pot.integrate_orbit(array(new_init), dt= -1, nsteps=nsteps)
	new_orbit = pot.integrate_orbit(reverse_orbit[-1], dt= 1, nsteps=nsteps*2)

	orbitg = new_orbit.to_coord_frame(coord.Galactic, galactocentric_frame=empty_galframe)
	orbit_orp = orbitg.transform_to(gc.Orphan)
	radvels = orbitg.radial_velocity

	orbit_orp.Lambda.wrap_angle=180.*u.deg

	f1 = interp1d(orbit_orp.Lambda, orbit_orp.Beta)
	f2 = interp1d(orbit_orp.Lambda, 5.*np.log10(orbit_orp.distance.value*1000.)-5.)
	f3 = interp1d(orbit_orp.Lambda, radvels)

	try: f1(ls)
	except ValueError: 
		#print 'f1 error'
		return -inf
	try: f2(ls)
	except ValueError: 
		#print 'f2 error'
		return -inf
	try: f3(ls)
	except ValueError: 
		#print 'f3 error'
		return -inf
	N_b = N(bs.value, f1(ls), (sb*np.ones(len(bs)))**2.)
	N_d = N(DMs, f2(ls), DMs_err**2. + sdm**2.)
	N_vr = N(vr.value, f3(ls), e_vr**2. + sv**2.)
	i_likelihood = N_b*N_d*N_vr

	if np.any(i_likelihood == 0):
		return -inf
	if np.any(i_likelihood == NaN):
		return -inf
	if np.any(i_likelihood == -NaN):
		return -inf
	#print np.sum(np.log(abs(i_likelihood)))
	return np.sum(np.log(abs(i_likelihood)))


def N(x, mu, sigma_sq):
	coeff = 1/np.sqrt(2*np.pi*sigma_sq)
	expo = - (x-mu)**2/(2*sigma_sq)
	return coeff * np.exp(expo)


def lnprior(p, pm_prior_on=pm_prior_on):
	b, dm, pm_l, pm_b, vrad, sb, sv, sdm, vc = p
	if sb <= 0:
		return -inf
	else: prior_sb = sb**-1
	if sv <= 0:
		return -inf
	else: prior_sv = sv**-1
	if sdm <= 0:
		return -inf
	else: prior_sdm = sdm**-1

	prior_dm = 10**(2*dm/5)

	if pm_prior_on:
		lpart = -(pm_l - (data_pm_l_cosb))**2. / (2*(0.05)**2)
		bpart = -(pm_b - (data_pm_b))**2. / (2*(0.05)**2)
		return lpart + bpart + np.log(prior_sb*prior_sv*prior_sdm*prior_dm)
	#else: 
	#	lpart = -(pm_l - (data_pm_l_cosb))**2. / (2*(0.25)**2)
	#	bpart = -(pm_b - (data_pm_b))**2. / (2*(0.25)**2)
	#	return lpart + bpart + np.log(prior_sb*prior_sv*prior_sdm*prior_dm)
	else: return np.log(prior_sb*prior_sv*prior_sdm*prior_dm)

def lnprob(p):
	lp = lnprior(p)
	if not np.isfinite(lp):
		return -inf
	return lp + log_likelihood(p)


p_start = [sohn_gal.b.value,5.*np.log10(newberg_interp_distance(sohn_gal.l.value)*1000.)-5.,sohn_gal.pm_l_cosb.value,sohn_gal.pm_b.value,float(newberg_interp_vrad(sohn_gal.l.value)), 3., 20., 5., 100.]

Nwalker,Ndim = 144,9
p0 = [p_start+[.2,.05,.1,.1,1.,1.,5.,1.,20.]*random.randn(Ndim) for i in range(Nwalker)]


sampler = emcee.EnsembleSampler(Nwalker,Ndim,lnprob, threads=16)
pos,prob,state = sampler.run_mcmc(p0, 5000)

#sampler.reset()
#pos,prob,state = sampler.run_mcmc(pos, 500)

samples = np.save(outname, (sampler.chain, sampler.flatchain))







