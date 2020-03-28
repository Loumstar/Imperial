import math as m
import numpy as np

from matplotlib import pyplot as plot
from matplotlib import animation

"""
String Oscillation methods
"""

def get_cfl(c, dt, dx):
    """
    Calculates the CFL of a finite difference method.
    CFL = c * (dt / dx)
    """
    return c * dt / dx

def string_is_invalid(cfl, u0, v0, nx):
    """
    Method to check if a numerical approximation is not
    convergent or is going to fail due to bad list sizes.
    """
    return cfl > 1 or len(v0) != nx or len(u0) != nx

def next_amplitude_1d(u, i, p, cfl):
    """
    Explicit finite difference method for calculating the amplitude
    of a point i along a string (1D) given that p > 1:

    u(i,p) = (2u(i,p-1) - u(i,p-2))
        + cfl^2 (u(i-1,p-1) - 2(i,p-1) - u(i+1,p-1))
    """
    finite_diff_time = (2 * u[p - 1][i]) - u[p - 2][i]
    finite_diff_x = u[p - 1][i - 1] - (2 * u[p - 1][i]) + u[p - 1][i + 1]
    
    return finite_diff_time + (cfl ** 2) * finite_diff_x

def start_amplitude_1d(u, i, dt, v0):
    """
    Explicit method for calculating the amplitude of a point i
    along a string (1D) at p = 1:

    u(i,p) = u(i,p-1) + dt u'(i,p-1)
    """
    return u[0][i] + dt * v0[i]

def string_frame(u, v0, p, bc, cfl, nx, dt):
    """
    Method to return the string amplitudes at a given frame.
    
    If a gridpoint is in the list of boundary conditions, the value of that
    boundary is added to the list. Therefore all boundary conditions must be
    Dirchlet type.
    
    If p = 1, use start_amplitude_1d(), else use normal explicit finite
    difference method next_amplitude_1d().
    """
    frame = []
    
    for i in range(nx):
        if i in bc:
            frame.append(0) # all boundary values are zero
        elif p == 1:
            frame.append(start_amplitude_1d(u, i, dt, v0))
        else:
            frame.append(next_amplitude_1d(u, i, p, cfl))

    return frame

def get_string_oscillation(length_string, t_end, u0, v0, bc, c, nx, nt):
    """
    Method to return a list of frames of the string oscillating.
    As each frame is a set of y values, the method returns a list of lists.

    Before looping through the frames, the set up is checked to see if it is
    invalid. If it is invalid, a list with the initial positions is returned.
    """
    u = [u0]
    dt = t_end / nt
    dx = length_string / nx

    cfl = get_cfl(c, dt, dx)

    if string_is_invalid(cfl, u0, v0, nx):
        print("Invalid string conditions")
        return u

    for p in range(1, nt):
        u.append(string_frame(u, v0, p, bc, cfl, nx, dt))
    
    return u

def find_bounds(i, bc, nx):
    """
    Method to find which boundary conditions the grid point i falls between.
    The method returns a lower bound and upper bound as grid points, in a tuple.
    """
    for j in range(1, len(bc)):
        if i >= bc[j - 1] and i <= bc[j]:
            return bc[j - 1], bc[j]

def velocity_from_bow(a, c, length, bc, nx, bow_location=0):
    """
    Method to return a list containing the velocity of each gridpoint along the
    string. The method tries to model a bow string the string at a point.

    The speed is calculated by finding the ratio of the position of each
    gridpoint to the length from the nearest bound to the bow:
    
    ratio = (i - lb) / (bow - lb) if i is between lb and bow.
    ratio = (ub - i) / (ub - bow) if i is between bow and ub.

    This list of ratios is then converted to a sinusoid to make the velocities
    more wave-like.

    If someone is holding down on a fret, the velocities between the neck and
    fret are zero, so any grid points that are not between the same boundary
    conditions as the bow will have zero velocity.
    """
    # velocity = angular speed * amplitude
    # angular speed = 2 * pi * f = 2 * pi * (wave speed / wavelength)
    # wavelength = 2 * string length because it is 1st mode of vibration
    max_velocity = a * (2 * m.pi * c / (2 * length))
    bow_gp = int(nx * bow_location / length)

    velocities = []
    
    for i in range(nx):
        lb, ub = find_bounds(i, bc, nx)
        if lb <= bow_gp and bow_gp <= ub: # if bow falls within these bounds
            ratio = (ub - i) / float(ub - bow_gp) if i >= bow_gp \
               else (i - lb) / float(bow_gp - lb)
            velocities.append(max_velocity * m.sin((m.pi / 2) * ratio))
        else:
            velocities.append(0)

    return velocities


def get_string_propagation(tension, mass_per_length):
    """
    Method to calculate the wave speed in a string.
    """
    return m.sqrt(tension / mass_per_length)

"""
Sound Oscillation Methods
"""

def air_is_invalid(cfl, a0, v0, ns):
    """
    Method to check if a numerical approximation is not
    convergent or is going to fail due to bad list sizes.
    """
    return cfl > 1 \
        or len(a0) != ns or len(a0[0]) != ns \
        or len(v0) != ns or len(v0[0]) != ns

def next_amplitude_2d(a, i, j, p, cfl):
    """
    Explicit finite difference method for calculating the amplitude
    of a point i,j in air (2D) given that p > 1:

    u(i,p) = (2u(i,j,p-1) - u(i,j,p-2))
        + cfl^2 (u(i-1,j,p-1) - 2(i,j,p-1) - u(i+1,j,p-1))
        + cfl^2 (u(i,j-1,p-1) - 2(i,j,p-1) - u(i,j+1,p-1))
    """
    finite_diff_time = (2 * a[p - 1][j][i]) - a[p - 2][j][i]
    finite_diff_x = a[p - 1][j][i - 1] - (2 * a[p - 1][j][i]) + a[p - 1][j][i + 1]
    finite_diff_y = a[p - 1][j - 1][i] - (2 * a[p - 1][j][i]) + a[p - 1][j + 1][i]
    
    return finite_diff_time + (cfl ** 2) * (finite_diff_x + finite_diff_y)

def start_amplitude_2d(a, i, j, dt, v0):
    """
    Explicit method for calculating the amplitude
    of a point i,j in air (2D) at p = 1:

    u(i,j,p) = u(i,j,p-1) + dt u'(i,j,p-1)
    """
    return a[0][j][i] + dt * v0[j][i]

def is_near_string(i, j, string_coord_gp, string_nx):
    """
    Method to determine whether a point i,j in air is near
    enough to the string to 'transfer' its amplitude.

    Returns a Boolean.
    """
    return i >= string_coord_gp[0] and i < string_coord_gp[0] + string_nx \
       and j == string_coord_gp[1]

def string_amplitude(i, string_coord_gp, string_frame):
    """
    Method to return the amplitude of the string at a point
    i in the air.
    """
    return string_frame[i - string_coord_gp[0]]

def air_frame(a, p, cfl, string, string_coord_gp, string_nx, v0, bcx, bcy, ns, ds, dt):
    """
    Method to return the air amplitudes at a given frame.
    
    If a gridpoint is in the list of boundary conditions, that gridpoint's
    value will be zero. Therefore all boundary conditions must be Dirchlet type.
    
    If a point i,j is near to the string:
        violin_x_coord < x < violin_x_coord + string_length
        y = violin_y_coord

    Then the amplitude of the string at that x coord is transferred to i,j.

    If p = 1, use start_amplitude_2d(), else use normal explicit finite
    difference method next_amplitude_2d().
    """
    frame = []
    
    for j in range(ns):
        row = []
        
        for i in range(ns):
            
            if i in bcx or j in bcy:
                row.append(0) # all boundary conditions are zero, simplifies code a little
            elif is_near_string(i, j, string_coord_gp, string_nx):
                row.append(string_amplitude(i, string_coord_gp, string[p]))
            elif p == 1:
                row.append(start_amplitude_2d(a, i, j, dt, v0))
            else:
                row.append(next_amplitude_2d(a, i, j, p, cfl))
        
        frame.append(row)

    return frame

def get_air_oscillation(length, t_end, a0, v0, string, string_coord, string_nx, c, bcx, bcy, ns, nt):
    """
    Method to return a list of frames of the air oscillating.
    As each frame is a set of x and y values, the method returns a list of list
    of lists.

    Before looping through the frames, the set up is checked to see if it is
    invalid. If it is invalid, a list with the initial positions is returned.
    """
    a = [a0]
    dt = t_end / nt
    ds = length / ns # ds is used to describe both dx and dy

    cfl = get_cfl(c, dt, ds)

    # Convert string coord to gridpoint indices
    string_x, string_y = string_coord
    string_x_gp = int(ns * string_x / length)
    string_y_gp = int(ns * string_y / length)

    string_coord_gp = (string_x_gp, string_y_gp)

    if air_is_invalid(cfl, a0, v0, ns):
        print("Invalid air conditions")
        return a

    # For each time step p
    for p in range(1, nt - 1):
        print("Air approximation: %.1f%% complete" % (100 * p / float(nt)))
        a.append(air_frame(a, p, cfl, string, string_coord_gp, string_nx, v0, bcx, bcy, ns, ds, dt))
    
    print("Air approximation: 100% complete")
    return a

def get_sound_propagation(gas_ratio, gas_constant, air_temperature):
    """
    Method to calculate the wave speed in air.
    """
    return m.sqrt(gas_ratio * gas_constant * (273.15 + air_temperature))

"""
Main Script.
Attempt to mimic open E on a violin.
"""

t_end = 0.008
nt = 800

dt = t_end / nt

"""
String properties and numerical approximation.
"""

string_nx = 50
string_length = 0.325

tension = 76.51
density = 7800
diameter = 0.26 * (10 ** -3)
amplitude = 5 * (10 ** -3)

mass_per_length = density * (diameter ** 2) * m.pi / 4
string_speed = get_string_propagation(tension, mass_per_length)

string_bc = [0, 49] # grid points at which the boundary value is 0

# String is flat, but has initial velocity from bow
string_u0 = np.zeros(string_nx)
string_v0 = velocity_from_bow(amplitude, string_speed, string_length, string_bc, string_nx, bow_location=0.300)

# Run approximation
string_x_values = np.linspace(0, string_length, num=string_nx, endpoint=True)
string_oscillation = get_string_oscillation(string_length, t_end, string_u0, string_v0,
                                                string_bc, string_speed, string_nx, nt)

"""
Air Properties and numerical approximaiton.
"""

room_ns = 400 # ns used to describe both nx and ny. Therefore room is a square
room_length = 2.6

speed_sound = get_sound_propagation(1.4, 287, 20)

# Position of leftmost grid point of string, relative to the top-left corner of the room.
string_coord = (1.13, 1.3) # Centre of the room

bcx_air = [0, 399] # x grid points at which the boundary value is 0
bcy_air = [0, 399] # y grid points at which the boundary value is 0

# Air is initial flat and still
air_u0 = np.zeros((room_ns, room_ns))
air_v0 = np.zeros((room_ns, room_ns))

# Run approximation
air_oscillation = get_air_oscillation(room_length, t_end, air_u0, air_v0, string_oscillation, string_coord,
                                                string_nx, speed_sound, bcx_air, bcy_air, room_ns, nt)

"""
Animate Plot.

Three animations:
    - Violin oscillation
    - Sound propagation through air
    - Sound detected at point (0, 0.65)
"""

print("Opening Matplotlib.")

fig1 = plot.figure()
fig1.canvas.set_window_title("Violin String")

fig2 = plot.figure()
fig2.canvas.set_window_title("Air Pressure")

fig3 = plot.figure()
fig3.canvas.set_window_title("Sound at (0, 0.65)")

ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
ax3 = fig3.add_subplot(111)

x_values = np.linspace(0, room_length, num=room_ns, endpoint=True)
x_axis, y_axis = np.meshgrid(x_values, x_values)

violin_images = []
air_images = []

sound_x_values = []
sound_y_values = []

sound_plot, = ax3.plot(sound_x_values, sound_y_values, 'k')

def init_sound():
    ax3.set_xlim(0, dt * len(air_oscillation))
    ax3.set_ylim(-0.005, 0.005)
    
    return sound_plot,

def update_sound(frame):
    sound_x_values.append(dt * frame)
    sound_y_values.append(air_oscillation[frame][99][199])
    sound_plot.set_data(sound_x_values, sound_y_values)
    
    return sound_plot,

for i, frame in enumerate(air_oscillation):
    violin_image = ax1.plot(string_x_values, string_oscillation[i], 'k-', animated=True)
    air_image = ax2.imshow(frame, vmin=-0.005, vmax=0.005, extent=[-1.3, 1.3, -1.3, 1.3],
                            aspect='auto', animated=True)

    violin_images.append(violin_image)
    air_images.append([air_image])

ani1 = animation.ArtistAnimation(fig1, violin_images, interval=50, blit=True)
ani2 = animation.ArtistAnimation(fig2, air_images, interval=50, blit=True)

ani3 = animation.FuncAnimation(fig3, update_sound, frames=range(len(air_oscillation)), init_func=init_sound, blit=True)

WriterClass = animation.writers['ffmpeg']
writer = WriterClass(fps=15, metadata=dict(artist='Louis Manestar'), bitrate=1800)

print("Writing Animations.")

ani1.save('violin_string.mp4', writer=writer)
ani2.save('violin_sound_waves.mp4', writer=writer)
ani3.save('violin_point_sound.mp4', writer=writer)

print("Done.")