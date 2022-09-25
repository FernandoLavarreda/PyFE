#Fernand José Lavarreda Urizar

import numpy as np
from typing import Tuple, List
from math import atan, pi, cos, sin

class Node:
    code = 0
    def __init__(self, x:float, y:float, name:str, limits:Tuple[bool, bool]=(False, False)):
        self.name = name
        self.x = x
        self.y = y
        self.restrictions = limits
        self.code = Node.code
        Node.code+=1
    
    
    def get_code(self):
        return self.code
    
    
    def get_name(self):
        return self.name
    
    
    def get_restrictions(self):
        return self.restrictions
    
    
    def get_x(self):
        return self.x
    
    
    def get_y(self):
        return self.y



class Element:
    code = 0
    def __init__(self, young_modulus:float, area:float, nodes:Tuple[Node, Node]):
        self.young = young_modulus
        self.area = area
        self.nodes = nodes
        self.angle = None
        self.length = None
        self.stiff_matrix = None
        self.trans_matrix = None
        self.thetha()
        self.length_()
        self.stiffness()
        self.transform()
        self.code = Element.code
        Element.code+=1
    
    
    def thetha(self):
        xcomp = self.nodes[1].get_x()-self.nodes[0].get_x()
        ycomp = self.nodes[1].get_y()-self.nodes[0].get_y()
        
        if xcomp < 0:
            self.angle = atan(ycomp/xcomp)+pi
        elif xcomp == 0:
            if ycomp < 0:
                self.angle = -pi/2
            else:
                self.angle = pi/2
        else:
            self.angle = atan(ycomp/xcomp)
    
    
    def stiffness(self):
        k = self.area*self.young/self.length
        matrix = np.array([[cos(self.angle)**2, sin(self.angle)*cos(self.angle), -cos(self.angle)**2, -sin(self.angle)*cos(self.angle)],\
                  [sin(self.angle)*cos(self.angle), sin(self.angle)**2, -sin(self.angle)*cos(self.angle), -sin(self.angle)**2],\
                  [-cos(self.angle)**2, -sin(self.angle)*cos(self.angle), cos(self.angle)**2, sin(self.angle)*cos(self.angle)],\
                  [-sin(self.angle)*cos(self.angle), -sin(self.angle)**2, sin(self.angle)*cos(self.angle), sin(self.angle)**2]])
        self.stiff_matrix = k*matrix
    
    
    def transform(self):
        self.trans_matrix = np.array([[cos(self.angle), sin(self.angle), 0, 0], [-sin(self.angle), cos(self.angle), 0, 0], [0, 0, cos(self.angle), sin(self.angle)], [0, 0, -sin(self.angle), cos(self.angle)]])
    
    def length_(self):
        xcomp = self.nodes[1].get_x()-self.nodes[0].get_x()
        ycomp = self.nodes[1].get_y()-self.nodes[0].get_y()
        
        self.length = ((xcomp**2+ycomp**2)**0.5)
    
    
    def get_transform(self):
        return self.trans_matrix
    
    
    def get_stiffness(self):
        return self.stiff_matrix
    
    
    def get_nodes(self):
        return self.nodes
    
    
    def get_length(self):
        return self.length
    
    
    def get_young(self):
        return self.young
    
    

def global_matrix(elements:List[Element]):
    gmatrix = np.zeros((Node.code*2, Node.code*2))
    for element in elements:
        codei = element.get_nodes()[0].get_code()
        codej = element.get_nodes()[1].get_code()
        for row in range(4):
            for column in range(4):
                if row < 2:
                    if column < 2:
                        gmatrix[codei*2+(row%2)][codei*2+(column%2)] += element.get_stiffness()[row][column]
                    else:
                        gmatrix[codei*2+(row%2)][codej*2+(column%2)] += element.get_stiffness()[row][column]
                else:
                    if column < 2:
                        gmatrix[codej*2+(row%2)][codei*2+(column%2)] += element.get_stiffness()[row][column]
                    else:
                        gmatrix[codej*2+(row%2)][codej*2+(column%2)] += element.get_stiffness()[row][column]
    return gmatrix


def solve_matrix(matrix:np.array, conditions:List[float], nodes:List[Node]):
    assert len(nodes) == len(matrix)//2
    boundary_cond = matrix.copy()
    counter = 0
    
    for node in nodes:
        if node.get_restrictions()[0]:
            line = [0 for i in range(len(nodes)*2)]
            line[counter*2] = 1
            boundary_cond[counter*2] = line
        
        if node.get_restrictions()[1]:
            line = [0 for i in range(len(nodes)*2)]
            line[counter*2+1] = 1
            boundary_cond[counter*2+1] = line
        counter+=1 
    
    inverse = np.linalg.inv(boundary_cond)
    displacement = np.matmul(inverse, conditions)
    reactions = np.matmul(matrix, displacement) - np.array(conditions)
    
    return displacement, reactions


def stress(elements:List[Element], displacement:np.array):
    local_disp = []
    stresses = []
    for element in elements:
        nodes = element.get_nodes()
        inode = nodes[0].get_code()
        jnode = nodes[1].get_code()
        global_disp = np.array([displacement[inode*2], displacement[inode*2+1], displacement[jnode*2], displacement[jnode*2+1]])
        locald = np.matmul(element.get_transform(), global_disp)
        local_disp.append(locald)
        stress = element.get_young()/element.get_length()*(locald[0]-locald[2])
        stresses.append(stress)
    return local_disp, stresses



class Node3D:
    code = 0
    def __init__(self, x, y, z, name, limits:Tuple[bool, bool, bool]=(False, False, False)):
        self.x = x
        self.y = y
        self.z = z
        self.name = name
        self.restrictions = limits
        self.code = Node3D.code
        Node3D.code+=1
    
    
    def get_name(self):
        return self.name
    
    
    def get_restrictions(self):
        return self.restrictions
    
    
    def get_code(self):
        return self.code
    
    
    def get_x(self):
        return self.x
    
    
    def get_y(self):
        return self.y
    
    
    def get_z(self):
        return self.z


class Element3D:
    code = 0
    def __init__(self, young_modulus:float, area:float, nodes:Tuple[Node3D, Node3D]):
        self.young = young_modulus
        self.area = area
        self.nodes = nodes
        self.length = None
        self.stiff_matrix = None
        self.length_()
        self.stiffness()
        self.code = Element3D.code
        Element3D.code+=1
    
    
    def length_(self):
        dx = (self.nodes[1].get_x()-self.nodes[0].get_x())**2
        dy = (self.nodes[1].get_y()-self.nodes[0].get_y())**2
        dz = (self.nodes[1].get_z()-self.nodes[0].get_z())**2
        self.length = (dx+dy+dz)**0.5
    
    
    def stiffness(self):
        k = self.area*self.young/self.length
        cx = (self.nodes[1].get_x()-self.nodes[0].get_x())/self.length
        cy = (self.nodes[1].get_y()-self.nodes[0].get_y())/self.length
        cz = (self.nodes[1].get_z()-self.nodes[0].get_z())/self.length
        matrix = np.array([[cx**2, cx*cy, cx*cz, -cx**2, -cx*cy, -cx*cz], [cx*cy, cy**2, cy*cz, -cx*cy, -cy**2, -cy*cz],\
                           [cx*cz, cz*cy, cz**2, -cx*cz, -cz*cy, -cz**2], [-cx**2, -cx*cy, -cx*cz, cx**2, cx*cy, cx*cz],\
                           [-cx*cy, -cy**2, -cy*cz, cx*cy, cy**2, cy*cz], [-cx*cz, -cz*cy, -cz**2, cx*cz, cz*cy, cz**2]])
        
        self.stiff_matrix = k*matrix
    
    def get_stiffness(self):
        return self.stiff_matrix
    
    def get_nodes(self):
        return self.nodes


def global_matrix3D(elements:List[Element3D]):
    gmatrix = np.zeros((Node3D.code*3, Node3D.code*3))
    for element in elements:
        codei = element.get_nodes()[0].get_code()
        codej = element.get_nodes()[1].get_code()
        for row in range(6):
            for column in range(6):
                if row < 3:
                    if column < 3:
                        gmatrix[codei*3+(row%3)][codei*3+(column%3)] += element.get_stiffness()[row][column]
                    else:
                        gmatrix[codei*3+(row%3)][codej*3+(column%3)] += element.get_stiffness()[row][column]
                else:
                    if column < 3:
                        gmatrix[codej*3+(row%3)][codei*3+(column%3)] += element.get_stiffness()[row][column]
                    else:
                        gmatrix[codej*3+(row%3)][codej*3+(column%3)] += element.get_stiffness()[row][column]
    return gmatrix


def solve_matrix3D(matrix:np.array, conditions:List[float], nodes:List[Node3D]):
    assert len(nodes) == len(matrix)//3
    boundary_cond = matrix.copy()
    counter = 0
    
    for node in nodes:
        if node.get_restrictions()[0]:
            line = [0 for i in range(len(nodes)*3)]
            line[counter*3] = 1
            boundary_cond[counter*3] = line
        
        if node.get_restrictions()[1]:
            line = [0 for i in range(len(nodes)*3)]
            line[counter*3+1] = 1
            boundary_cond[counter*3+1] = line
        
        if node.get_restrictions()[2]:
            line = [0 for i in range(len(nodes)*3)]
            line[counter*3+2] = 1
            boundary_cond[counter*3+2] = line
        
        counter+=1 
    
    inverse = np.linalg.inv(boundary_cond)
    displacement = np.matmul(inverse, conditions)
    reactions = np.matmul(matrix, displacement) - np.array(conditions)
    
    return displacement, reactions


if __name__ == "__main__":
    """
    node1 = Node(0, 0, "i", (True, True))
    node2 = Node(36, 0, "j")
    node3 = Node(0, 36, "j", (True, True))
    node4 = Node(36, 36, "j")
    node5 = Node(72, 36, "j")
    
    sad1 = Element(1.9e6, 8, (node1, node2))
    sad2 = Element(1.9e6, 8, (node2, node3))
    sad3 = Element(1.9e6, 8, (node3, node4))
    sad4 = Element(1.9e6, 8, (node2, node4))
    sad5 = Element(1.9e6, 8, (node2, node5))
    sad6 = Element(1.9e6, 8, (node4, node5))
    #print(sad1.get_stiffness())
    nodes = [node1, node2, node3, node4, node5]
    elss = [sad1, sad2, sad3, sad4, sad5, sad6]
    mat = global_matrix(elss)
    r, r1 = solve_matrix(mat, [0, 0, 0, 0, 0, 0, 0, -500, 0, -500], nodes)
    localsd, strs = stress(elss, r)
    print(strs)
    #print(r)
    
    #########3
    #Prueba 3d
    nodos = [Node3D(0, 0, 0, "1"), Node3D(-2, 0, -1.5, "2", (True, True, True)), Node3D(-2, 0, 1.5, "3", (True, True, True)), Node3D(-2, 1.5, 0, "4", (True, True, True))] 
    elementos = [Element3D(70e9, 15/100**2, (nodos[0], nodos[2])), Element3D(70e9, 15/100**2, (nodos[0], nodos[1])), Element3D(70e9, 15/100**2, (nodos[0], nodos[3]))]
    
    matrizg = global_matrix3D(elementos)
    global_desplazamiento, reacciones = solve_matrix3D(matrizg, [0, -5000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], nodos)
    print(reacciones)
    print(global_desplazamiento)
    ##########33
    """
    
    #Problema de clase para armaduras
    nodos = [Node(0, 0, "1", (True, True)), Node(96, 0, "2"), Node(192, 0, "3", (False, True)), Node(48, 6.9282*12, "4"), Node(144, 6.9282*12, "5")]
    elementos = [Element(30e6, 3, (nodos[0], nodos[1])), Element(30e6, 3, (nodos[1], nodos[2])), Element(30e6, 3, (nodos[0], nodos[3])), \
                 Element(30e6, 3, (nodos[3], nodos[4])), Element(30e6, 3, (nodos[4], nodos[2])),\
                 Element(30e6, 3, (nodos[3], nodos[1])), Element(30e6, 3, (nodos[1], nodos[4]))]
    
    matrizg = global_matrix(elementos)
    global_desplazamiento, reacciones = solve_matrix(matrizg, [0, 0, 0, -5000, 0, 0, 0, 0, 0, 0], nodos)
    locales_desplazamiento, esfuerzos = stress(elementos, global_desplazamiento)
    
    counter = 0
    for esfuerzo in esfuerzos:
        print(f"Esfuerzo en elemento {counter+1} que conecta los nodos({elementos[counter].get_nodes()[0].get_code()+1}, {elementos[counter].get_nodes()[1].get_code()+1}) es: {round(-esfuerzo, 2)}psi")
        counter+=1
    
    for i in range(1, len(nodos)*2+1):
        if i %2 == 0:
            r = "y"
        else:
            r = "x"
        print("R"+r+f" para nodo ({(i+1)//2}) es: {round(reacciones[i-1])}lb")
        print(f"Desplazamiento {r} del nodo ({(i+1)//2}) es: {round(global_desplazamiento[i-1], 6)}in")
    print("\n")
    
    
    #Problema en 3D
    nodos = [Node3D(48, 0, 0, "A"), Node3D(0, 0, -24, "B", (True, True, True)), Node3D(0, 0, 24, "C", (True, True, True)), Node3D(0, 48, 0, "D", (True, True, True)), Node3D(0, 0, 0, "E")]
    elementos = [Element3D(29e6, 3.093, (nodos[0], nodos[1])), Element3D(29e6, 3.093, (nodos[0], nodos[4])), Element3D(29e6, 3.093, (nodos[0], nodos[2])),\
                 Element3D(29e6, 3.093, (nodos[0], nodos[3])), Element3D(29e6, 3.093, (nodos[1], nodos[3])), Element3D(29e6, 3.093, (nodos[2], nodos[3])),\
                 Element3D(29e6, 3.093, (nodos[4], nodos[3])), Element3D(29e6, 3.093, (nodos[4], nodos[2])), Element3D(29e6, 3.093, (nodos[4], nodos[1]))]
    
    matrizg = global_matrix3D(elementos)
    global_desplazamiento, reacciones = solve_matrix3D(matrizg, [0, -250, -150, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -250, -150], nodos)
    pos = ["x", "y", "z"]
    for i in range(len(nodos)*3):
        r = pos[i%3]
        print(f"Desplazamiento del nodo {nodos[i//3].get_name()} en {r} es: {round(global_desplazamiento[i], 7)}in")
    
    
    
    
    
    
    
    
    
    
    