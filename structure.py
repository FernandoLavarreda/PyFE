from typing import List
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from elementos import Node, Node3D, Element, Element3D


def structure2D(elements:List[Element], nodes:List[Node]):
    figure = plt.figure()
    axes = figure.add_subplot(111)
    for node in nodes:
        axes.scatter(node.get_x(), node.get_y(), color="r")
    for element in elements:
        node1, node2 = element.get_nodes()
        axes.plot([node1.get_x(), node2.get_x()], [node1.get_y(), node2.get_y()], color="black")
    return figure


def structure3D(elements:List[Element], nodes:List[Node3D]):
    figure = plt.figure()
    axes = Axes3D(figure)
    figure.add_axes(axes)
    for node in nodes:
        axes.scatter3D(node.get_x(), node.get_y(), node.get_z(), color="red")
    for element in elements:
        node1, node2 = element.get_nodes()
        axes.plot3D([node1.get_x(), node2.get_x()], [node1.get_y(), node2.get_y()], [node1.get_z(), node2.get_z()], color="black")
    return figure


if __name__ == "__main__":
    nodos = [Node(0, 0, "1", (True, True)), Node(96, 0, "2"), Node(192, 0, "3", (False, True)), Node(48, 6.9282*12, "4"), Node(144, 6.9282*12, "5")]
    elementos = [Element(30e6, 3, (nodos[0], nodos[1])), Element(30e6, 3, (nodos[1], nodos[2])), Element(30e6, 3, (nodos[0], nodos[3])), \
                 Element(30e6, 3, (nodos[3], nodos[4])), Element(30e6, 3, (nodos[4], nodos[2])),\
                 Element(30e6, 3, (nodos[3], nodos[1])), Element(30e6, 3, (nodos[1], nodos[4]))]
    
    structure2D(elementos, nodos)
    plt.show()
    
    nodos = [Node3D(48, 0, 0, "A"), Node3D(0, 0, -24, "B", (True, True, True)), Node3D(0, 0, 24, "C", (True, True, True)), Node3D(0, 48, 0, "D", (True, True, True)), Node3D(0, 0, 0, "E")]
    elementos = [Element3D(29e6, 3.093, (nodos[0], nodos[1])), Element3D(29e6, 3.093, (nodos[0], nodos[4])), Element3D(29e6, 3.093, (nodos[0], nodos[2])),\
                 Element3D(29e6, 3.093, (nodos[0], nodos[3])), Element3D(29e6, 3.093, (nodos[1], nodos[3])), Element3D(29e6, 3.093, (nodos[2], nodos[3])),\
                 Element3D(29e6, 3.093, (nodos[4], nodos[3])), Element3D(29e6, 3.093, (nodos[4], nodos[2])), Element3D(29e6, 3.093, (nodos[4], nodos[1]))]
    
    structure3D(elementos, nodos)
    plt.show()

