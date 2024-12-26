import shapely
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.plotting import plot_polygon, plot_points

points = shapely.MultiPoint([(0.0, 0.0), (1.0, 1.0)])
s = points.buffer(2.0)

s = shapely.Polygon(s.buffer(1.0).exterior)
print(s)
fig = plt.figure(1,  dpi=90)
plot_polygon(s)
plt.show()

#
# # 2: complex line
# fig = plt.figure(1,  dpi=90)
# line2 = LineString([(0, 0), (1, 1), (0, 2), (2, 2), (-1, 1), (1, 0)])
#
#
# plot_line(line2, ax=ax, add_points=False, color=YELLOW, alpha=0.7)
# plot_points(line2, ax=ax, color=GRAY, alpha=0.7)
# plot_points(line2.boundary, ax=ax, color=BLACK)
#
# ax.set_title('b) complex')
#
# set_limits(ax, -2, 3, -1, 3)
#
# plt.show()
