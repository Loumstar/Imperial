import math as m
import numpy as np

from violin_class import ViolinOscillation

class StringSoundPropagation:
    def __init__(self, x_length, y_length, **kwargs):
        self.dx = 0
        self.dy = 0
        self.dt = 0
        
        self.cfl = m.inf

        self.__x_length = 1.0
        self.__y_length = 1.0
        self.__t_end = 1.0

        self.__nx = 1
        self.__ny = 1
        self.__nt = 1

        self.__c = 0.0

        self.x_length = x_length
        self.y_length = y_length

        self.nx = kwargs.get('nx', 1)
        self.ny = kwargs.get('ny', 1)
        self.nt = kwargs.get('nt', 1)

        self.c = kwargs.get('c', 0.0)

        self.u = []

        self.string = kwargs.get('string', ViolinOscillation(1.0))
        self.string_coord = kwargs.get('string_coord', (0.0, 0.0))

        self.u0 = kwargs.get('u0', np.zeros((self.ny, self.nx)))
        self.v0 = kwargs.get('v0', np.zeros((self.ny, self.nx)))

    """
    Setters and Getters for Simulation Properties.
    If one changes, the CFL and other properties are recalculated.
    """

    @property
    def nx(self):
        return self.__nx

    @nx.setter
    def nx(self, nx):
        if type(nx) != int:
            print("Number of gridpoints must be an integer.")
        else:
            self.__nx = nx
            self.dx = self.x_length / self.nx
            self.cfl = self.calculate_cfl()

    @property
    def ny(self):
        return self.__ny

    @ny.setter
    def ny(self, ny):
        if type(ny) != int:
            print("Number of gridpoints must be an integer.")
        else:
            self.__ny = ny
            self.dy = self.y_length / self.ny
            self.cfl = self.calculate_cfl()

    @property
    def x_length(self):
        return self.__x_length

    @x_length.setter
    def x_length(self, x_length):
        if x_length == 0:
            print("Length cannot be zero.")
        else:
            self.__x_length = float(x_length)
            self.dx = self.x_length / self.nx
            self.cfl = self.calculate_cfl()

    @property
    def y_length(self):
        return self.__y_length

    @y_length.setter
    def y_length(self, y_length):
        if y_length == 0:
            print("Length cannot be zero.")
        else:
            self.__y_length = float(y_length)
            self.dy = self.y_length / self.ny
            self.cfl = self.calculate_cfl()

    @property
    def nt(self):
        return self.__nt

    @nt.setter
    def nt(self, nt):
        if type(nt) != int:
            print("Number of samples must be an integer.")
        else:
            self.__nt = nt
            self.dt = self.t_end / self.nt
            self.cfl = self.calculate_cfl()

    @property
    def t_end(self):
        return self.__t_end

    @t_end.setter
    def t_end(self, t_end):
        if t_end == 0:
            print("Warning: simulation time set to zero.")
        
        self.__t_end = float(t_end)
        self.dt = self.t_end / self.nt
        self.cfl = self.calculate_cfl()

    @property
    def c(self):
        return self.__c

    @c.setter
    def c(self, c):
        self.__c = float(c)
        self.cfl = self.calculate_cfl()

    @property
    def u0(self):
        return self.__u0

    @u0.setter
    def u0(self, u0):
        if len(u0) != self.ny or len(u0[0]) != self.nx:
            print("Initial velocity length must equal the number of gridpoints.")
        else:
            self.__u0 = u0

    @property
    def v0(self):
        return self.__v0

    @v0.setter
    def v0(self, v0):
        if len(v0) != self.ny or len(v0[0]) != self.nx:
            print("Initial velocity length must equal the number of gridpoints.")
        else:
            self.__v0 = v0

    @property
    def string(self):
        return self.__string

    @string.setter
    def string(self, string):
        if not isinstance(string, ViolinOscillation):
            print("Must be a ViolinOscillation instance.")
        else:
            self.__string = string

    @property
    def string_coord(self):
        return self.__string_coord

    @string_coord.setter
    def string_coord(self, string_coord):
        if string_coord[0] >= 0 and string_coord[0] < self.x_length and \
           string_coord[1] >= 0 and string_coord[1] < self.y_length:
            self.__string_coord = string_coord
        else:
            print("String falls outside simulation boundaries.")   

    """
    Methods for Simulation Properties.
    """

    def __amplitude(self, i, j, p):
        fd_t = (2 * self.u[p - 1][j][i]) - self.u[p - 2][j][i]
        fd_x = self.u[p - 1][j][i - 1] - (2 * self.u[p - 1][j][i]) + self.u[p - 1][j][i + 1]
        fd_y = self.u[p - 1][j - 1][i] - (2 * self.u[p - 1][j][i]) + self.u[p - 1][j + 1][i]
        
        return fd_t + ((self.cfl ** 2) * (fd_x + fd_y))

    def __second_frame_amplitude(self, i, j):
        return self.u[0][j][i] + (self.dt * self.v0[j][i])

    def __nearest_string_value(self, i, j, p):
        xgp = int(((self.dx * i) - self.string_coord[0]) / self.string.dx)
        
        if xgp > 0 and xgp < self.string.nx:
            y = (self.dy * j) - self.string_coord[1] + self.string.u[p][xgp]
            
            if y > 0 and y < self.dy:
                return True, self.string.u[p][xgp]
        
        return False, 0.0

    def __pressure_frame(self, p):
        frame = []

        for j in range(self.ny):
            row = []
            
            for i in range(self.nx):
                (near_string, string_amplitude) = self.__nearest_string_value(i, j, p)

                if i == 0 or i == self.nx - 1 \
                or j == 0 or j == self.ny - 1:
                    row.append(0)
                elif p == 1:
                    row.append(self.__second_frame_amplitude(i, j))
                elif near_string:
                    row.append(string_amplitude)
                else:
                    row.append(self.__amplitude(i, j, p))
            
            frame.append(row)

        return frame

    def __pressure_distribution(self):
        self.u.append(self.u0)
        
        for p in range(1, self.nt):
            print("Sound simulation %.1f%% complete." % (100 * p / float(self.nt)))
            self.u.append(self.__pressure_frame(p))

        print("Sound simulation 100.0% complete.")

    """
    User methods for simulation.
    """

    def reset_data(self):
        self.u = []

    def set_propogation_speed(self, gas_ratio, gas_constant, air_temperature):
        self.c = m.sqrt(gas_ratio * gas_constant * (air_temperature + 273.15))

    def calculate_cfl(self):
        if self.dx == 0 or self.dy == 0:
            return m.inf

        cfl_x = self.c * self.dt / self.dx
        cfl_y = self.c * self.dt / self.dy

        return cfl_x if cfl_x > cfl_y else cfl_y

    def is_invalid(self):
        return self.cfl > 1 \
            or self.string.dt != self.dt \
            or len(self.string.u) != self.nt \
            or len(self.u0[0]) != self.nx \
            or len(self.u0) != self.ny \
            or len(self.v0[0]) != self.nx \
            or len(self.v0) != self.ny

    def get_pressure_distribution(self, t_end):
        self.reset_data()
        self.t_end = t_end

        if self.is_invalid():
            print("Invalid simulation conditions")
            print("CFL = %.1f" % self.cfl)
            return self.u

        self.__pressure_distribution()
        
        return self.u

    def get_pressure_at_point(self, x, y):
        pressure = []
        
        if len(self.u) == 0:
            print("No simulation data.")
            return pressure
        elif x < 0 or x > self.x_length or y < 0 or y > self.y_length:
            print("Point coordinates falls outside simulation boundaries")
            return pressure
        
        i = int(x / self.dx)
        j = int(y / self.dy)

        for frame in self.u:
            pressure.append(frame[j][i])

        return pressure