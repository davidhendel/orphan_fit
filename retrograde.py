import cld_fnc as cf
from astropy.io import ascii
from astropy.coordinates import ICRS, Galactic, SkyCoord, Distance
from astropy import units as u
import astropy.coordinates as coord
import gary.coordinates as gc
from gary.units import galactic
from scipy.interpolate import interp1d
import scipy
import gary
import astropy.constants as ac
from scipy.interpolate import interp1d
import emcee

if 0:
	fdakldjalsk=True
	#v_gsr=v_helio+10.1*cos(b)*cos(l)+ 224*cos(b)*sin(l)+6.7*sin(b)

	#Format Units    Label   Explanations
	#-----------------------------------------------------------------------
	#A4     ---      Name    RR Lyrae identifier
	#F10.6  deg      RAdeg   Right Ascension in decimal degrees (J2000)
	#F9.6   deg      DEdeg   Declination in decimal degrees (J2000)
	#A6     ---      Survey  Survey
	#F4.1   kpc      HDis    Heliocentric distance
	#F6.3   mag      <mag>   Flux-averaged magnitude (1)
	#F5.3   mag      Ar      SDSS r band extinction (2)
	#F4.2   mag      Amp     Amplitude
	#F5.2   mag      mag0    Magnitude at epoch of maximum brightness
	#F8.6   d        Per     Period of pulsation
	#F12.6  d        HJD     Reduced Heliocentric Julian Date of maximum
	#                         brightness; HJD-2400000
	#F6.1   km/s     VHel    Heliocentric velocity
	#F4.1   km/s   e_VHel    Uncertainty in VHel
	#F6.1   km/s     VGSR    Galactic standard of rest velocity; Section 3
	#F5.2   mas/yr   pml     Proper motion along Galactic longitude
	#F5.2   mas/yr   pmb     Proper motion along Galactic latitude
	#F4.2   mas/yr   epm     Error in each proper motion
	#F5.2   [-]      [Fe/H]  Metallicity
	#A6     ---      Mem     Probability of being a Orphan stream member
	#-----------------------------------------------------------------------

#make Newberg potential
pot = gary.potential.CompositePotential()
pot['bulge'] = gary.potential.HernquistPotential(m=3.4e10, c=0.7, units=(u.kpc,u.Myr,u.Msun,u.rad))
pot['disk']  = gary.potential.MiyamotoNagaiPotential(m=1.0e11,a=6.5,b=0.26, units=(u.kpc,u.Myr,u.Msun,u.rad))
pot['halo']  = gary.potential.LogarithmicPotential(v_c=(73.*np.sqrt(2.)*u.km/u.second).to(u.kpc/u.Myr).value, r_h=12.,q1=1., q2=1., q3=1., units=(u.kpc,u.Myr,u.Msun,u.rad))

#load data
ls, bs, dbs, vgsrs, dvgsrs, distss, ddistss = np.loadtxt('/Users/hendel/projects/orphan/newberg_data.txt', usecols=(1,2,3,6,7,10,11), skiprows=2, unpack=True)

rrly_datafile = '/Users/hendel/projects/orphan/sesar_rrly.txt'
rrly_data = ascii.read(rrly_datafile)
rrly_data = rrly_data[rrly_data["Mem"]=='high']

#RRLy coords
rr_c = SkyCoord(ra=rrly_data['RAdeg'], dec = rrly_data['DEdeg'], distance=rrly_data["HDis"], unit=(u.degree,u.degree,u.kpc))
rr_cg = rr_c.galactic
#F turnoff cords
f_cg = Galactic(b=bs, l=ls, unit=(u.degree, u.degree), distance=Distance(distss, u.kpc))
f_vgsrs = vgsrs

#newberg ICs
(l,b,R)=(218, 53.5, 28.6)
newberg_ic = Galactic(l = l, b= b, unit=(u.degree, u.degree), distance=Distance(R,u.kpc))

newberg_init = [newberg_ic.cartesian.x.value-8., 
		newberg_ic.cartesian.y.value, 
		newberg_ic.cartesian.z.value, 
		(-156/1000.), 
		(79/1000.), 
		(107/1000.)]

reverse_orbit = pot.integrate_orbit(array(newberg_init), dt= -.1, nsteps=11000)
orbit = pot.integrate_orbit(np.append(reverse_orbit.pos[:,-1].value, reverse_orbit.vel[:,-1].value), dt= .1, nsteps=22000)

orbitg, vels = orbit.to_frame(coord.Galactic)
radvels = vels[2].to(u.km/u.s)
gsrvels = gc.vhel_to_vgsr(orbitg, radvels)


grillmair_c = SkyCoord(ra=167, dec = -14, distance=18., unit=(u.degree,u.degree, u.kpc))
grillmair_cg = grillmair_c.galactic
grillmair_mEncl = pot.mass_enclosed([grillmair_cg.cartesian.x.value-8,grillmair_cg.cartesian.y.value,grillmair_cg.cartesian.z.value])
galcen_dist = np.sqrt((grillmair_cg.cartesian.x.value-8)**2.+grillmair_cg.cartesian.y.value**2. + grillmair_cg.cartesian.z.value**2.)

obs_pos = SkyCoord(ra=167*np.ones(1000), dec = -14*np.ones(1000), distance=np.arange(13,23,.01), unit=(u.degree,u.degree, u.kpc))
obs_posg=obs_pos.galactic
obs_xyz = obs_posg.cartesian


#find nearest orbit point

orb_dist = np.sqrt(
	(orbit.pos[0,:].value-(grillmair_cg.cartesian.x.value-8.))**2 + 
	(orbit.pos[1,:].value-grillmair_cg.cartesian.y.value)**2    + 
	(orbit.pos[2,:].value-grillmair_cg.cartesian.z.value)**2    )
closest = where(orb_dist == min(orb_dist))[0]


L_x = (orbit.pos[1,closest]*orbit.vel[2,closest] - orbit.pos[2,closest]*orbit.vel[1,closest])
L_y = (orbit.pos[2,closest]*orbit.vel[0,closest] - orbit.pos[0,closest]*orbit.vel[2,closest])
L_z = (orbit.pos[0,closest]*orbit.vel[1,closest] - orbit.pos[1,closest]*orbit.vel[0,closest])
L_mag_orbit = np.sqrt(L_x**2. + L_y**2. + L_z**2.)

E_orbit = 0.5*(orbit.vel[0,closest]**2 + orbit.vel[0,closest]**2 + orbit.vel[0,closest]**2) + pot.value([orbit.pos[0,closest],orbit.pos[1,closest],orbit.pos[2,closest]])*u.kpc*u.kpc/u.Myr/u.Myr

r_tide = 1*u.kpc
r_p = 15*u.kpc
e_s = 2*r_tide*grillmair_mEncl*u.Msun*ac.G/(galcen_dist**2*u.kpc**2)
e_s=e_s.to(u.kpc*u.kpc/u.Myr/u.Myr)

l_s = (sqrt(3)+2)*L_mag_orbit*(r_tide/r_p)

L_mag_sat = L_mag_orbit-l_s*.5
E_sat = E_orbit-e_s*2.

v_sat = np.sqrt(2*E_sat.value - pot.value([(grillmair_cg.cartesian.x.value-8)**2.+grillmair_cg.cartesian.y.value**2. + grillmair_cg.cartesian.z.value**2.]))*u.kpc/u.Myr
satvs = orbit.vel[:,closest]*v_sat.value/norm(orbit.vel[:,closest].value)

grillmair_init= [grillmair_cg.cartesian.x.value-8.,
	grillmair_cg.cartesian.y.value,
	grillmair_cg.cartesian.z.value,
	satvs[0].value,
	satvs[1].value,
	satvs[2].value]

reverse_orbit_sat = pot.integrate_orbit(array(grillmair_init), dt= -.1, nsteps=11000)
orbit_sat = pot.integrate_orbit(np.append(reverse_orbit_sat.pos[:,-1].value, reverse_orbit_sat.vel[:,-1].value), dt= .1, nsteps=22000)

plt.plot(np.sqrt(orbit_sat.pos[0,:]**2+orbit_sat.pos[1,:]**2+orbit_sat.pos[2,:]**2))

#load sim stuff
if 1:
	# sim stuff
	#2 apo ago, 3 2.98 Gyr
	# <Quantity [ 55.12632742,-27.96102663,-54.53747397] kpc>
	# <Quantity [ 0.00553293,-0.06312098,-0.00862088] kpc / Myr>

	#3 apo ago, 4.9573 Gyr
	#<Quantity [-18.49800083,-69.503532  , -6.4507776 ] kpc>
	#<Quantity [-0.06344375,-0.05863569, 0.04447514] kpc / Myr>

	d = cf.particles()
	d.init('/Users/hendel/projects/orphan/nbody/2k3/', 40)

	xyz = coord.Galactocentric([d.x, d.y, d.z] * u.kpc, galcen_distance=8.*u.kpc)
	vxyz = [d.vx, d.vy, d.vz] * u.km/u.s
	d.xyzg = xyz.transform_to(coord.Galactic)
	pm_ra,pm_dec,vr = gc.vgal_to_hel(d.xyzg, vxyz)
	d.nbodygsrvels = gc.vhel_to_vgsr(d.xyzg, vr)


vgsr_diffs = np.zeros(len(rr_cg.l.value))
metals = np.zeros(len(rr_cg.l.value))
for i in np.arange(len(rr_cg.l.value)):
	l_pos = where(abs(orbitg.l.value-rr_cg.l.value[i]) == np.min(abs(rr_cg.l.value[i] - orbitg.l.value)))

	vgsr_diffs[i] = abs(gsrvels.value[l_pos] -rrly_data['VGSR'].data[i])/rrly_data['e_VHel'].data[i]
	metals[i]= rrly_data["[Fe/H]"].data[i]






[-0.85193276 -0.05140298]


#############
newberg_init = [newberg_ic.cartesian.x.value-8., 
		newberg_ic.cartesian.y.value, 
		newberg_ic.cartesian.z.value, 
		(-156/1000.), 
		(79/1000.), 
		(107/1000.)]

c = coord.Galactocentric(newberg_ic.cartesian.x - 8*u.kpc, 
		newberg_ic.cartesian.y, 
		newberg_ic.cartesian.z)
newberg_vxyz = 		[-156., 79., 107.]*u.km/u.s
pm_l, pm_b, vrad = gc.vgal_to_hel(c.transform_to(coord.Galactic), newberg_vxyz)
pm_l_norm = pm_l.value/(np.sqrt(pm_l.value**2+pm_b.value**2))
pm_b_norm = pm_b.value/(np.sqrt(pm_l.value**2+pm_b.value**2))
norm = (np.sqrt(pm_l.value**2+pm_b.value**2))

theta= 0.
[pm_l_norm_rot, pm_b_norm_rot] = np.dot([pm_l_norm, pm_b_norm], [[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

# p = mag, theta

def orb_var(mag, theta, nsteps = 400):

	[pm_l_norm_rot, pm_b_norm_rot] = np.dot([pm_l_norm, pm_b_norm], [[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
	new_pms = [pm_l_norm_rot*mag, pm_b_norm_rot*mag]*u.mas/u.yr
	new_xyz_vels =  gc.vhel_to_gal(c.transform_to(coord.Galactic), pm = new_pms, rv = vrad)

	new_init = [newberg_ic.cartesian.x.value-8., 
			newberg_ic.cartesian.y.value, 
			newberg_ic.cartesian.z.value, 
			(new_xyz_vels[0].value/1000.), 
			(new_xyz_vels[1].value/1000.), 
			(new_xyz_vels[2].value/1000.)]

	reverse_orbit = pot.integrate_orbit(array(new_init), dt= -1, nsteps=nsteps)
	new_orbit = pot.integrate_orbit(np.append(reverse_orbit.pos[:,-1].value, reverse_orbit.vel[:,-1].value), dt= 1, nsteps=nsteps*2)

	orbitg, vels = orbit.to_frame(coord.Galactic)
	radvels = vels[2].to(u.km/u.s)
	gsrvels = gc.vhel_to_vgsr(orbitg, radvels)

	return new_orbit, orbitg, gsrvels


def log_likelihood(p, nsteps=500, ls = rr_cg.l, bs= rr_cg.b, dists = rr_cg.distance.value, 
				vgsr = rrly_data["VHel"], e_decs = 1., e_dists= rr_cg.distance*.1, e_vgsr = rrly_data["e_VHel"]):

	mag, theta = p
	[pm_l_norm_rot, pm_b_norm_rot] = np.dot([pm_l_norm, pm_b_norm], [[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
	new_pms = [pm_l_norm_rot*mag, pm_b_norm_rot*mag]*u.mas/u.yr
	new_xyz_vels =  gc.vhel_to_gal(c.transform_to(coord.Galactic), pm = new_pms, rv = vrad)

	new_init = [newberg_ic.cartesian.x.value-8., 
			newberg_ic.cartesian.y.value, 
			newberg_ic.cartesian.z.value, 
			(new_xyz_vels[0].value/1000.), 
			(new_xyz_vels[1].value/1000.), 
			(new_xyz_vels[2].value/1000.)]

	reverse_orbit = pot.integrate_orbit(array(new_init), dt= -1, nsteps=nsteps)
	new_orbit = pot.integrate_orbit(np.append(reverse_orbit.pos[:,-1].value, reverse_orbit.vel[:,-1].value), dt= 1, nsteps=nsteps*2)

	orbitg, vels = new_orbit.to_frame(coord.Galactic)
	radvels = vels[2].to(u.km/u.s)
	gsrvels = gc.vhel_to_vgsr(orbitg, radvels)

	f1 = interp1d(orbitg.l, orbitg.b)
	f2 = interp1d(orbitg.l, orbitg.distance)
	f3 = interp1d(orbitg.l, gsrvels.to(u.km/u.second))
	try: f1(ls)
	except ValueError: return -inf
	N_b = N(bs.value, f1(ls), 3*np.ones(len(bs)))
	N_d = N(dists, f2(ls), dists*.1)
	N_vgsr = N(vgsr.data, f3(ls), e_vgsr.data)
	i_likelihood = N_b*N_d*N_vgsr

	print p
	print np.sum(np.log(abs(i_likelihood)))
	return np.sum(np.log(abs(i_likelihood)))


def N(x, mu, sigma):
	coeff = 1/np.sqrt(2*np.pi*sigma**2)
	expo = - (x-mu)**2/(2*sigma**2)
	return coeff * np.exp(expo)


def lnprior(p):
	mag, theta = p
	if abs(mag) >3.:
		return -inf
	return 0

def lnprob(p):
	lp = lnprior(p)
	if not isfinite(lp):
		return -inf
	return lp + log_likelihood(p)


p_start = [norm, 0.]

Nwalker,Ndim = 30,2
p0 = [p_start+.1*random.randn(Ndim) for i in range(Nwalker)]


sampler = emcee.EnsembleSampler(Nwalker,Ndim,lnprob, threads=3)
pos,prob,state = sampler.run_mcmc(p0, 500)

sampler.reset()
pos,prob,state = sampler.run_mcmc(pos, 500)


m_mag, m_theta = median(sampler.flatchain, axis=0)

new_orbit, new_orbitg, new_vgsrs = orb_var(m_mag, m_theta)
xyz_plot(new_orbit)
plt.figure()
bvr_plot(new_orbit)
plt.figure()
plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
corner.corner(sampler.flatchain, labels=['mag', 'theta'])

################

for i in np.arange(0,30):

	factor = np.arange(-3, 0, .1)*norm
	new_pms = [pm_l_norm_rot*factor[i], pm_b_norm_rot*factor[i]]*u.mas/u.yr
	new_xyz_vels =  gc.vhel_to_gal(c.transform_to(coord.Galactic), pm = new_pms, rv = vrad)

	new_init = [newberg_ic.cartesian.x.value-8., 
			newberg_ic.cartesian.y.value, 
			newberg_ic.cartesian.z.value, 
			(new_xyz_vels[0].value/1000.), 
			(new_xyz_vels[1].value/1000.), 
			(new_xyz_vels[2].value/1000.)]

	reverse_orbit = pot.integrate_orbit(array(new_init), dt= -1, nsteps=1100)
	new_orbit = pot.integrate_orbit(np.append(reverse_orbit.pos[:,-1].value, reverse_orbit.vel[:,-1].value), dt= 1, nsteps=2200)
	xyz_plot(new_orbit, color = cm.jet(i*8), data=False, sim=None)
	print i

xyz_plot(orbit)
xyz_plot(orbit_sat, color='cyan')

################


def xyz_plot(orbit, sim=None, data=True, color='k', newfig=False):
	if newfig==True: 
		plt.figure()

	plt.subplot(1,3,1, aspect='equal')
	if sim != None: 
		plt.scatter(sim.x, sim.y, alpha=0.25, c='k', s=2, edgecolor='none')
	if data==True:
		#Sesar RR Lyrae 
		plt.scatter(rr_cg.cartesian.x.value-8., rr_cg.cartesian.y.value, c='g', alpha = 0.7, edgecolor='none')
		#Newberg BHB 
		plt.scatter(f_cg.cartesian.x.value-8., f_cg.cartesian.y.value, c='r', alpha = 0.7, edgecolor='none')
		#Grillmair projenitor
		plt.scatter(grillmair_cg.cartesian.x.value-8., grillmair_cg.cartesian.y.value, c='r', marker='*', s=150, zorder=10, edgecolor='none')
		#Grillmair position distance uncertianty
		plt.plot(obs_xyz.x.value-8., obs_xyz.y.value, c='r', lw=1.5)


	plt.plot(orbit.pos[0,:], orbit.pos[1,:],lw=2, c=color)

	plt.xlim([-70,70])
	plt.ylim([-40,100])
	plt.xlabel('x')
	plt.ylabel('y')


	plt.subplot(1,3,2, aspect='equal')
	if sim != None: 
		plt.scatter(sim.z, sim.y, alpha=0.25, c='k', s=2, edgecolor='none')
	if data==True:
		#Sesar RR Lyrae 
		plt.scatter(rr_cg.cartesian.z.value, rr_cg.cartesian.y.value, c='g', alpha = 0.7, edgecolor='none')
		#Newberg BHB 
		plt.scatter(f_cg.cartesian.z.value, f_cg.cartesian.y.value, c='r', alpha = 0.7, edgecolor='none')
		#Grillmair projenitor
		plt.scatter(grillmair_cg.cartesian.z.value, grillmair_cg.cartesian.y.value, c='r', marker='*', s=150, zorder=10, edgecolor='none')
		#Grillmair position distance uncertianty
		plt.plot(obs_xyz.z.value, obs_xyz.y.value, c='r', lw=1.5)


	plt.plot(orbit.pos[2,:], orbit.pos[1,:],lw=2, c=color)

	plt.xlim([-70,70])
	plt.ylim([-40,100])
	plt.xlabel('z')
	plt.ylabel('y')

	plt.subplot(1,3,3, aspect='equal')
	if sim != None: 
		plt.scatter(sim.x, sim.z, alpha=0.25, c='k', s=2, edgecolor='none')
	if data==True:
		#Sesar RR Lyrae 
		plt.scatter(rr_cg.cartesian.x.value-8., rr_cg.cartesian.z.value, c='g', alpha = 0.7, edgecolor='none')
		#Newberg BHB 
		plt.scatter(f_cg.cartesian.x.value-8., f_cg.cartesian.z.value, c='r', alpha = 0.7, edgecolor='none')
		#Grillmair projenitor
		plt.scatter(grillmair_cg.cartesian.x.value-8., grillmair_cg.cartesian.z.value, c='r', marker='*', s=150, zorder=10, edgecolor='none')
		#Grillmair position distance uncertianty
		plt.plot(obs_xyz.x.value-8., obs_xyz.z.value, c='r', lw=1.5)


	plt.plot(orbit.pos[0,:], orbit.pos[2,:],lw=2, c=color)

	plt.xlim([-70,70])
	plt.ylim([-70,70])
	plt.xlabel('x')
	plt.ylabel('z')


def bvr_plot(orbit, sim=None, data=True, color='k', newfig=False):
	
	orbitg, vels = orbit.to_frame(coord.Galactic)
	radvels = vels[2].to(u.km/u.s)
	gsrvels = gc.vhel_to_vgsr(orbitg, radvels)


	if newfig==True: 
		plt.figure()

	plt.subplot(3,1,1)
	if sim != None:
		plt.scatter(sim.xyzg.l.value, sim.xyzg.b.value, color='k', s=3, alpha=0.3, edgecolor='none')

	if data ==True:
		#Sesar RR Lyrae
		plt.scatter(rr_cg.l.value, rr_cg.b.value, c ='g', s = 15, alpha = 1, edgecolor='k', zorder=8)
		#Newberg BHB 
		plt.scatter(f_cg.l.value, f_cg.b.value, c='r',s=15, alpha = 1, edgecolor='k', zorder =9)
		#Grillmair projenitor
		plt.scatter(grillmair_cg.l.value, grillmair_cg.b.value, c='r', marker='*', s=150, zorder=10, edgecolor='k')
		#Grillmair position distance uncertianty
		plt.plot(obs_posg.l.value, obs_posg.b.value, c='r', lw=1.5)

	plt.plot(orbitg.l.value, orbitg.b.value, c = color, linewidth=3)

	plt.xlim([300,120])
	plt.ylim([0,70])
	plt.gca().set_xticklabels('',visible=False)
	plt.ylabel('b [degrees]')

	plt.subplot(3,1,2)
	if sim != None:
		plt.scatter(sim.xyzg.l.value, sim.nbodygsrvels.value, color='k', s=3, alpha=0.3, edgecolor='none')

	if data ==True:
		#Sesar RR Lyrae
		plt.scatter(rr_cg.l.value, rrly_data['VGSR'].data, c ='g', s = 15, alpha = 1, edgecolor='k', zorder=8)
		#Newberg BHB 
		plt.scatter(f_cg.l.value, f_vgsrs, c='r', alpha = 1, edgecolor='k', zorder=9)
		#Grillmair projenitor
		#unknown 
			#plt.scatter(grillmair_cg.l.value, grillmair_cg.b.value, c='r', marker='*', s=150, zorder=10, edgecolor='none')
		#Grillmair position distance uncertianty
		#unknown 
			#plt.plot(obs_posg.l.value, obs_posg.b.value, c='r', lw=1.5)
	
	plt.plot(orbitg.l.value, gsrvels, c = color, linewidth=3)

	plt.xlim([300,120])
	plt.ylim([-100,200])
	plt.gca().set_xticklabels('',visible=False)
	plt.ylabel('$v_{gsr}\ [km\ s^{-1}]$')


	plt.subplot(3,1,3)

	if sim != None:
		plt.scatter(sim.xyzg.l.value, sim.xyzg.distance.value, color='k', s=3, alpha=0.3, edgecolor='none')

	if data ==True:
		#Sesar RR Lyrae
		plt.scatter(rr_cg.l.value, rr_cg.distance.value, c ='g', s = 15, alpha = 1, edgecolor='k', zorder =8)
		#Newberg BHB 
		plt.scatter(f_cg.l.value, f_cg.distance.value, c='r', alpha = 1, edgecolor='k', zorder=9)
		#Grillmair projenitor
		plt.scatter(grillmair_cg.l.value, grillmair_cg.distance.value, c='r', marker='*', s=150, zorder=10, edgecolor='k')
		#Grillmair position distance uncertianty
		plt.plot(obs_posg.l.value, obs_posg.distance.value, c='r', lw=1.5)

	#plt.errorbar(ls,distss, yerr=ddistss, c='r', marker ='o', ls="None")
	sel = [(orbitg.l.value > 120) & (orbitg.l.value < 300)]
	plt.plot(orbitg[sel].l.value, orbitg[sel].distance.value, c = color, linewidth=3)

	plt.xlim([300,120])
	plt.ylim([0,70])
	plt.xlabel('l [degrees]')
	plt.ylabel('Heliocentric distance [kpc]')



