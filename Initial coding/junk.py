# import platform
# print(platform.system())

# import re

# a = "A"
# b = ".*"
# m = re.search(b,a )
# print(m)
# import array as myarray

# first = myarray.array('b',[20,30,40])
# print(first)

# first.remove(30)
# print(first)

# Program to display calendar of the given month and year

# # importing calendar module
# import calendar
# import datetime

# currentDateTime = datetime.datetime.now()
# date = currentDateTime.date()

# # yy = date.year  # year
# # mm = date.month    # month

# To take month and year input from the user
# yy = int(input("Enter year: "))
# mm = int(input("Enter month: "))

# # display the calendar
# print(calendar.month(yy, mm))

# num = int(input("Enter number: "))
# if num != 0 and num != 2 and num != 5 and num not in range(10, 16) and num not in range(20, 26):
#     print('Okay')

x0 = 500
y0 = 0
z0 = 0
u = 0
b = 0.0002
a = 0.0095
l = 0.0001
d = 0.0001
x = 0
y = 0
z = 0

x_evo = []
y_evo = []
z_evo = []
for i in range(0,500):
    x = x0 + u - a*x0*y0 - b*x0
    y = y0 + a*x0*y0 + l*z0 - d*x0*y0
    z = z0 + b*x0 + d*x0*y0 - l*z0
    x_evo.append(x)
    y_evo.append(y)
    z_evo.append(z)
    print("x: %s, y: %s, z:%s" % (x, y, z))
    x0 = x
    y0 = y
    z0 = z

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# create data
time = list(range(500))
# plot
plt.plot(time, x_evo)
plt.plot(time, y_evo)
plt.plot(time, z_evo)
plt.xlim(right=10)
plt.ylim(bottom=-1, top=2)
plt.show()

