from math import pi, cos, sin
from numpy import array,dot, hstack, eye, vstack, empty, arctan2, sqrt

# ПРЯМАЯ ЗАДАЧА


def ht(q,d,a,alpha):
    if q == pi/2 or q == -pi/2:
        c_q = 0
    else:
        c_q = cos(q)
    if alpha == pi/2 or alpha == -pi/2:
        c_alpha = 0
    else:
        c_alpha = cos(alpha)
    T1 = array([[c_q, -sin(q), 0, 0], [sin(q), c_q, 0, 0], [0, 0, 1, 0],[0, 0, 0, 1]])
    T2 = vstack((hstack((eye(3),array([[0],[0],[d]]))),array([[0, 0, 0, 1]])))
    T3 = vstack((hstack((eye(3), array([[a], [0], [0]]))), array([[0, 0, 0, 1]])))
    T4 = array([[1, 0, 0, 0], [0, c_alpha, -sin(alpha), 0], [0, sin(alpha), c_alpha, 0],[0, 0, 0, 1]])
    return dot(dot(dot(T1,T2),T3),T4)


# Параметризация матрицы поворотов с помощью углов Эйлера
def eul(R):
    phi, theta, psi = 0, 0, 0
    if abs(R[2,2]) < 1:
        phi = arctan2(R[1,2],R[0,2])
        theta = arctan2(sqrt(1 - R[2,2] ** 2), R[2,2])
        psi = arctan2(R[2,1], -R[2,0])
    elif R[2,2] == 1:
        phi = 0
        theta = 0
        psi = arctan2(R[1,0], R[0,0])
    elif R[2,2] == -1:
        phi = arctan2(-R[0,1], -R[0,0])
        theta = pi
        psi = 0
    Ang=[round(phi,4),round(theta,4),round(psi,4)]
    return Ang


# Задание параметров Денавита-Хартенберга
q = array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
q0 = array([0, 0, pi/2, 0, 0, 0])
q = q + q0
a = array([0, 1, 0, 0, 0, 0])
d = array([1, 0, 0, 1, 0, 1])
alpha = array([pi/2, 0, pi/2, -pi/2, pi/2, 0])

T = empty([len(q), 16]) # пустая матрица Т размером 6 на 16, все 6 матриц представить как строку
for i in range(len(q)):
    T[i, :] = ht(q[i], d[i], a[i], alpha[i]).reshape(-1) # расчёт и представление в виде строки

T02 = T[0,:].reshape(4,4).dot(T[1,:].reshape(4,4)) # умножение Т01 на Т12
T03 = T02.dot(T[2,:].reshape(4,4))
T04 = T03.dot(T[3,:].reshape(4,4))
T05 = T04.dot(T[4,:].reshape(4,4))
T06 = T05.dot(T[5,:].reshape(4,4))

R = T06[0:3,0:3]

Koord = [round(T06[0,3],4),round(T06[1,3],4),round(T06[2,3],4)]

xi = Koord + eul(R)

print(xi)


# ОБРАТНАЯ ЗАДАЧА

def R06(phi,theta,psi):
    R06 = dot([[cos(phi), -sin(phi), 0],[sin(phi), cos(phi), 0],[0, 0, 1]],[[cos(theta), 0, sin(theta)],[0, 1, 0],[-sin(theta), 0, cos(theta)]])
    R06 = dot(R06,[[cos(psi), -sin(psi), 0],[sin(psi), cos(psi), 0],[0, 0, 1]])
    return R06


def kord(R06,d,a):
    p06 = array([[x],[y],[z]])
    p04 = p06 - d[5] * dot(R06,[[0],[0],[1]])

    xc, yc, zc = [float(p04[i]) for i in range(3)]
    q[0]= round(arctan2(yc,xc),1)
    cosq3 = ((zc - d[0]) ** 2 + xc ** 2 + yc ** 2 - a[1] ** 2 - d[3] ** 2) / (2 * a[1] * d[3])

    if int(float(cosq3) / 2) == 1:
        q[1] = round((arctan2((zc - d[0]), sqrt(xc ** 2 + yc ** 2))),1)
        q[2] = 0
    elif int(float(cosq3) / 2) == -1:
        q[2] = pi
    elif abs(int(float(cosq3) / 2)) < 1:
        q[2] = round(arctan2(sqrt(1 - cosq3 ** 2), cosq3),1)
    q[1] = round(arctan2(zc - d[0], sqrt(xc ** 2 + yc ** 2)) - arctan2(d[3] * sin(q[2]), a[1] + d[3] * cos(q[2])),1)
    return q


def R(q,d,a,alpha,R06):
    T01 = ht(q[0],d[0],a[0], alpha[0])
    T12 = ht(q[1],d[1],a[1], alpha[1])
    T23 = ht(q[2] + pi / 2, d[0], a[2], alpha[2])
    T03 = dot(dot(T01,T12), T23)
    R03 = T03[0:3, 0:3]
    return dot(R03.T, R06)


# задание начальных условий
q = [float(0) for i in range(6)]
x, y, z, phi, theta, psi = [xi[i] for i in range(len(xi))]
a = array([0, 1, 0, 0, 0, 0])
d = array([1, 0, 0, 1, 0, 1])
alpha = array([pi/2, 0, pi/2, -pi/2, pi/2, 0])


qk = kord(R06(phi, theta, psi), d, a)
q = qk + eul(R(q, d, a, alpha, R06(phi, theta, psi)))
print(q)



