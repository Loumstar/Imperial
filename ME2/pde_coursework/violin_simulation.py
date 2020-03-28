import math as m
import numpy as np

from matplotlib import pyplot as plot
from matplotlib import animation

from violin_class import ViolinOscillation
from sound_class import StringSoundPropagation

tension = 76.51
density = 7800
string_length = 0.325
diameter = 0.26 * (10 ** -3)
amplitude = 5 * (10 ** -3)

string_nx = 50

gas_ratio = 1.4
gas_constant = 287.0
air_temp = 20.0

mass_per_length = density * (diameter ** 2) * m.pi / 4

violin = ViolinOscillation(string_length, nx=50, nt=800)
air = StringSoundPropagation(2.6, 2.6, nx=400, ny=400, nt=800, string=violin, string_coord=[1.15, 1.3])

violin.set_propogation_speed(tension, mass_per_length)
violin.set_initial_velocity(amplitude, bow_position=0.30)

air.set_propogation_speed(gas_ratio, gas_constant, air_temp)

oscillation = violin.get_string_oscillation(0.008)
pressure = air.get_pressure_distribution(0.008)
sound = air.get_pressure_at_point(1.3, 0.65)

print("Opening Matplotlib.")

figure1 = plot.figure()
figure1.canvas.set_window_title("Violin String")

figure2 = plot.figure()
figure2.canvas.set_window_title("Sound Waves")

figure3 = plot.figure()
figure3.canvas.set_window_title("Sound at (0, 0.65)")

ax1 = figure1.add_subplot(111)
ax2 = figure2.add_subplot(111)
ax3 = figure3.add_subplot(111)

violin_images = []
pressure_images = []

sound_x = []
sound_y = []

sound_plot, = ax3.plot(sound_x, sound_y, 'k')

def init_sound():
    ax3.set_xlim(0, air.t_end)
    ax3.set_ylim(-0.005, 0.005)
    return sound_plot,

def update_sound(p):
    sound_x.append(air.dt * p)
    sound_y.append(sound[p])
    sound_plot.set_data(sound_x, sound_y)
    
    return sound_plot,

x_values = np.linspace(0, string_length, num=50, endpoint=True)

for i, frame in enumerate(oscillation):
    violin_image = ax1.plot(x_values, frame, 'k', animated=True)
    pressure_image = ax2.imshow(pressure[i], extent=[-1.3, 1.3, -1.3, 1.3], vmin=-amplitude, vmax=amplitude, aspect='auto', animated=True)

    violin_images.append(violin_image)
    pressure_images.append([pressure_image])

ani1 = animation.ArtistAnimation(figure1, violin_images, interval=50)
ani2 = animation.ArtistAnimation(figure2, pressure_images, interval=50)

ani3 = animation.FuncAnimation(figure3, update_sound, frames=range(len(sound)), init_func=init_sound)

WriterClass = animation.writers['ffmpeg']
writer = WriterClass(fps=15, metadata=dict(artist='Louis Manestar'), bitrate=1800)

print("Writing Animations.")

ani1.save('violin_string.mp4', writer=writer)
ani2.save('violin_sound_waves.mp4', writer=writer)
ani3.save('violin_point_sound.mp4', writer=writer)

print("Done.")
