from scipy.optimize import minimize, rosen, rosen_der
import numpy as np 

# def func(x):
#     return -1/((x-3)**2+5)

# def der(x):
#     return 2*(x-3)/(((x-3)**2+5)**2)

# def func(x):
#     E= -1/((x[0]-x[1]-3)**2 +5)
    
#     return E
        


# def der(x):
#     x_der = [0,0]
#     x_der[0] = 2*(x[0]-x[1]-3)/(((x[0]-x[1]-3)**2+5)**2)
#     x_der[1] = -2*(x[0]-x[1]-3)/(((x[0]-x[1]-3)**2+5)**2)
#     print("x_der[0]: {}".format(x_der[0]))
#     print("x_der[1]: {}".format(x_der[1]))
#     return x_der


def func(x):
    E = 0
    for i in range(len(x)-1):
        for j in range(i+1, len(x)):
            if i < j:
                E += -1/((x[i]-x[j]-3)**2+5)
    print("E: {}".format(E))
    return E

def der(x):
    vec_der = list()
    for i in range(len(x)):
        sum_der = 0
        for j in range(len(x)):
            if i < j:
                sum_der += 2*(x[i]-x[j]-3)/(((x[i]-x[j]-3)**2+5)**2)
            elif i > j:
                sum_der += -2*(x[j]-x[i]-3)/(((x[j]-x[i]-3)**2+5)**2)
        vec_der.append(sum_der)
        print("{}: {}".format(i, sum_der))
    return vec_der


x= [10, 40]
res = minimize(func, x, method="BFGS", jac = der, options={'disp':True})
print("BFGS: {}".format(res.x))

