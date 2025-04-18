def find_angles(r1, r2, r3, r4, theta2):
    x = np.arange(-180, 180, 0.001) * np.pi / 180
    theta2_rad = theta2 * np.pi / 180

    f1 = np.arccos((r2 * np.cos(theta2_rad) + r3 * np.cos(x) - r1) / r4)
    f2 = -f1
    g1 = np.arcsin((r2 * np.sin(theta2_rad) + r3 * np.sin(x)) / r4)
    g2 = np.pi - g1

    diff1 = np.abs(f1 - g1)
    diff2 = np.abs(f1 - g2)
    diff3 = np.abs(f2 - g1)
    diff4 = np.abs(f2 - g2)

    tolerance = 0.0001
    indices1 = np.where(diff1 < tolerance)[0]
    indices2 = np.where(diff2 < tolerance)[0]
    indices3 = np.where(diff3 < tolerance)[0]
    indices4 = np.where(diff4 < tolerance)[0]

    intersections_x = np.concatenate((x[indices1], x[indices2], x[indices3], x[indices4]))
    intersections_y = np.concatenate((f1[indices1], f1[indices2], f2[indices3], f2[indices4]))

    plt.plot(x, f1, label='f1(x)')
    plt.plot(x, g1, label='g1(x)')
    plt.plot(x, f2, label='f2(x)')
    plt.plot(x, g2, label='g2(x)')
    plt.plot(intersections_x, intersections_y, 'ro', label='Intersections')
    plt.legend()
    plt.xlabel('x (radians)')
    plt.ylabel('angle (radians)')
    plt.title('Intersections of f(x) and g(x)')
    plt.grid()
    plt.show()

    index_of_separation = len(intersections_x)
    for i in range(1, len(intersections_x)):
        if intersections_x[i] - intersections_x[i - 1] > 0.1:
            index_of_separation = i
            break

    theta3 = intersections_x[(index_of_separation + len(intersections_x) - 1) // 2]
    theta4 = intersections_y[(index_of_separation + len(intersections_x) - 1) // 2]

    return theta3, theta4


def find_angular_velocity(r1, r2, r3, r4, theta2, omega2, theta3, theta4):
    theta2_rad = theta2 * np.pi / 180

    omega4 = ((r2 * omega2 * np.cos(theta2_rad) - (1 / np.tan(theta3)) * r2 * omega2 * np.sin(theta2_rad)) /
              (r4 * np.cos(theta4) - r4 * np.sin(theta4) * (1 / np.tan(theta3))))
    omega3 = (r4 * omega4 * np.sin(theta4) - r2 * omega2 * np.sin(theta2_rad)) / (r3 * np.sin(theta3))

    return omega3, omega4

def find_angular_acceleration(r1, r2, r3, r4, theta2, omega2, alpha2, theta3, theta4, omega3, omega4):
    theta2_rad = theta2 * np.pi / 180

    x = symbols('x')

    a = (r3 * omega3 ** 2 * np.cos(theta3) + r3 * x * np.sin(theta3) +
         r2 * omega2 ** 2 * np.cos(theta2_rad) + r2 * alpha2 * np.sin(theta2_rad) -
         r4 * np.cos(theta4) * omega4 ** 2) / (r4 * np.sin(theta4))

    b = (-r3 * omega3 ** 2 * np.sin(theta3) + r3 * x * np.cos(theta3) -
         r2 * omega2 ** 2 * np.sin(theta2_rad) + r2 * alpha2 * np.cos(theta2_rad) +
         r4 * np.sin(theta4) * omega4 ** 2) / (r4 * np.cos(theta4))

    equation = Eq(a, b)
    solution = solve(equation, x)

    alpha3 = solution[0]
    alpha4 = a.subs(x, alpha3)

    return alpha3, alpha4