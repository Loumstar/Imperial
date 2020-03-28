import math as m

class ViolinOscillation:
    def __init__(self, length, **kwargs):
        self.dx = 0
        self.dt = 0

        self.cfl = m.inf

        self.__length = 1.0
        self.__t_end = 1.0

        self.__nx = 2
        self.__nt = 2

        self.__c = 0.0

        self.length = length

        self.nx = kwargs.get('nx', 2)
        self.nt = kwargs.get('nt', 2)

        self.c = kwargs.get('c', 0.0)

        self.u = []

        self.bc = kwargs.get('bc', [])

        self.u0 = kwargs.get('u0', [0.0] * self.nx)
        self.v0 = kwargs.get('v0', [0.0] * self.nx)

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
            self.dx = self.length / self.nx
            self.cfl = self.calculate_cfl()

    @property
    def length(self):
        return self.__length

    @length.setter
    def length(self, length):
        if length == 0:
            print("Length cannot be zero.")
        else:
            self.__length = float(length)
            self.dx = self.length / self.nx
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
    def bc(self):
        return self.__bc

    @bc.setter
    def bc(self, bc):
        if type(bc) != list:
            print("Boundary conditions must be a list.")
        elif len(bc) > 0 and max(bc) >= self.nx:
            print("Boundary points must be in the range (0, %i)." % (self.nx - 1))
        else:
            self.__bc = bc

    @property
    def u0(self):
        return self.__u0

    @u0.setter
    def u0(self, u0):
        if len(u0) != self.nx:
            print("Initial velocity length must equal the number of gridpoints.")
        else:
            self.__u0 = u0

    @property
    def v0(self):
        return self.__v0

    @v0.setter
    def v0(self, v0):
        if len(v0) != self.nx:
            print("Initial velocity length must equal the number of gridpoints.")
        else:
            self.__v0 = v0

    """
    Methods for Simulation Properties.
    """

    def __find_bounds(self, i):
        if len(self.bc) == 0:
            return 0, self.nx - 1
        elif i <= self.bc[0]:
            return 0, self.bc[0]
        elif i > self.bc[-1]:
            return self.bc[-1], self.nx - 1

        for j in range(1, len(self.bc)):
            if i > self.bc[j - 1] and i <= self.bc[j]:
                return self.bc[j - 1], self.bc[j]

    def __velocity_from_bow(self, a, bow_gp):
        max_velocity = a * (m.pi * self.c / self.length)
        velocities = []
        
        for i in range(self.nx):
            lb, ub = self.__find_bounds(i)
            if lb <= bow_gp and bow_gp < ub:
                ratio = (ub - i) / float(ub - bow_gp) if i >= bow_gp \
                   else (i - lb) / float(bow_gp - lb)
                velocities.append(max_velocity * m.sin((m.pi / 2) * ratio))
            else:
                velocities.append(0)

        return velocities

    def __amplitude(self, i, p):
        fd_t = (2 * self.u[p - 1][i]) - self.u[p - 2][i]
        fd_x = self.u[p - 1][i - 1] - (2 * self.u[p - 1][i]) + self.u[p - 1][i + 1]
        
        return fd_t + ((self.cfl ** 2) * fd_x)

    def __second_amplitude(self, i):
        return self.u[0][i] + (self.dt * self.v0[i])

    def __string_frame(self, p):
        frame = []

        for i in range(self.nx):
            if i == 0 or i == self.nx - 1:
                frame.append(0)
            elif i in self.bc:
                frame.append(0)
            elif p == 1:
                frame.append(self.__second_amplitude(i))
            else:
                frame.append(self.__amplitude(i, p))

        return frame

    def __string_oscillation(self):
        self.u.append(self.u0)
        
        for p in range(1, self.nt):
            print("Violin simulation %.1f%% complete." % (100 * p / float(self.nt)))
            self.u.append(self.__string_frame(p))

        print("Violin simulation 100.0% complete.")

    """
    User methods for simulation.
    """

    def reset_data(self):
        self.u = []

    def set_propogation_speed(self, tension, mass_per_length):
        self.c = m.sqrt(tension / mass_per_length)

    def calculate_cfl(self):
        if self.dx == 0:
            return m.inf
        else:
            return self.c * self.dt / self.dx

    def is_invalid(self):
        return self.cfl > 1 \
            or len(self.v0) != self.nx \
            or len(self.u0) != self.nx

    def set_initial_velocity(self, a, bow_position=0):
        bow_gp = int(self.nx * bow_position / self.length) - 1
        self.v0 = self.__velocity_from_bow(a, bow_gp)

    def get_string_oscillation(self, t_end):
        self.reset_data()
        self.t_end = t_end

        if self.is_invalid():
            print("Invalid simulation conditions")
            print("CFL = %.1f" % self.cfl)
            return self.u

        self.__string_oscillation()
        
        return self.u
