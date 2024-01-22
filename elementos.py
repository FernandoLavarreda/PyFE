#!./venv/bin/python
#Fernand José Lavarreda Urizar

import numpy as np
from typing import Tuple, List
from dataclasses import dataclass
from math import atan, pi, cos, sin

class Node:
    """2D Vector representation has a name and movement in X and/or Y direction may be restricted
    defaulted to false in both directions"""
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
    """Model composed of two 2D Nodes, has uniform cross section as Young-modulus"""
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
        """Determine rotation of Element respect global coordinate system"""
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
    """Element matrix to solve"""
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
    """Obtain displacement of nodes and global reactions"""
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
    """Obtain stresses from global displacement"""
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
    """3D Vector representation, has a name and movement in X, Y and/or Z direction may be restricted
    defaulted to false in all directions"""
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
    """Model made of 2 3D Nodes"""
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
    remove_cond = []
    conditions_copy = conditions.copy()
    for node in nodes:
        if node.get_restrictions()[0]:
            line = [0 for i in range(len(nodes)*3)]
            line[counter*3] = 1
            boundary_cond[counter*3] = line
            remove_cond.append(counter*3)
            
        if node.get_restrictions()[1]:
            line = [0 for i in range(len(nodes)*3)]
            line[counter*3+1] = 1
            boundary_cond[counter*3+1] = line
            remove_cond.append(counter*3+1)
            
        if node.get_restrictions()[2]:
            line = [0 for i in range(len(nodes)*3)]
            line[counter*3+2] = 1
            boundary_cond[counter*3+2] = line
            remove_cond.append(counter*3+2)
            
        counter+=1 
    
    for i in remove_cond:
        conditions[i] = 0
    
    inverse = np.linalg.inv(boundary_cond)
    displacement = np.matmul(inverse, conditions)
    reactions = np.matmul(matrix, displacement) - np.array(conditions_copy)
    return displacement, reactions



class Node3DoF(Node):
    """Reinterpretation of Node3D, instead of using the 3rd dimension solely as Z movement, this model
    may be used to interpret Frames in 2D where the elements may rotate"""
    def __init__(self, x:float, y:float, name:str, limits:Tuple[bool, bool, bool]=(False, False, False)):
        super().__init__(x, y, name, limits)


class ElementFrame(Element):
    def __init__(self, young_modulus:float, area:float, inertia:float, nodes:Tuple[Node, Node]):
        super().__init__(young_modulus, area, nodes)
        self.inertia = inertia
        self.trans_transpose = None
        self.transformf()
        self.stiffnessf()
        
    
    
    def stiffnessf(self):
        k = self.area*self.young/self.length
        sub = 12*self.young*self.inertia/self.length**3
        sub2 = 6*self.young*self.inertia/self.length**2 
        sub3 = 4*self.young*self.inertia/self.length 
        sub4 = 2*self.young*self.inertia/self.length 
        
        local_stiffness = np.array([[k, 0, 0, -k, 0, 0], [0, sub, sub2, 0, -sub, sub2], [0, sub2, sub3, 0, -sub2, sub4],\
                                    [-k, 0, 0, k, 0, 0], [0, -sub, -sub2, 0, sub, -sub2], [0, sub2, sub4, 0, -sub2, sub3]])
        self.stiff_matrix = self.trans_transpose@local_stiffness@self.trans_matrix
    
    
    def transformf(self):
        self.trans_matrix = np.array([[cos(self.angle), sin(self.angle), 0, 0, 0, 0], [-sin(self.angle), cos(self.angle), 0, 0, 0, 0],\
                                      [0, 0, 1, 0, 0, 0], [0, 0, 0, cos(self.angle), sin(self.angle), 0],\
                                      [0, 0, 0, -sin(self.angle), cos(self.angle), 0], [0, 0, 0, 0, 0, 1]])
        self.trans_transpose = np.transpose(self.trans_matrix)


def globalize_perpendicular_load(load:float, angle:float):
        """Obtain x and y component of force perpendicular to given angle"""
        xcomp = cos(angle+pi/2)*load
        ycomp = sin(angle+pi/2)*load
        return xcomp, ycomp


@dataclass
class PerpendicularLoad:
    magnitude:float = 0
    reactioni:Tuple[float, float] = (0 , 0)
    reactionj:Tuple[float, float] = (0 , 0)
    momenti:float = 0
    momentj:float = 0
    
    def apply(self, element):
        return
    

class PointLoad(PerpendicularLoad):
    def __init__(self, magnitude, distance):
        """distance refers to the distance from node i of the element
        """
        assert distance >= 0, "Distance must be positive"
        self.magnitude = magnitude
        self.distance = distance
    
    
    def apply(self, element):
        assert self.distance <= element.length, "Cannot apply load outside element"
        reactioni = -self.magnitude*((element.length-self.distance)/element.length)**2*(3-2*(element.length-self.distance)/element.length)
        reactionj = -self.magnitude*(self.distance/element.length)**2*(3-2*self.distance/element.length)
        self.momenti = -self.magnitude*self.distance*((element.length-self.distance)/element.length)**2
        self.momentj = self.magnitude*(self.distance/element.length)**2*(element.length-self.distance)
        self.reactioni = globalize_perpendicular_load(reactioni, element.angle)
        self.reactionj = globalize_perpendicular_load(reactionj, element.angle)
    

class UniformLoad(PerpendicularLoad):
    def __init__(self, magnitude):
        self.magnitude = magnitude
    
    
    def apply(self, element):
        reactioni = -self.magnitude*element.length/2
        reactionj = reactioni
        self.momenti = -self.magnitude*element.length**2/12
        self.momentj = -self.momenti
        self.reactioni = globalize_perpendicular_load(reactioni, element.angle)
        self.reactionj = globalize_perpendicular_load(reactionj, element.angle)


class TriangleLoad(PerpendicularLoad):
    def __init__(self, magnitude, ascending=True):
        """
        Ascending indicates positive slope (increase in absolute value of the force)
        """
        self.magnitude = magnitude
        self.ascending = ascending
    
    def apply(self, element):
        reactioni = -3*element.length*self.magnitude/20
        reactionj = -7*element.length*self.magnitude/20
        self.momenti = -self.magnitude*element.length**2/30
        self.momentj = self.magnitude*element.length**2/20
        
        if not self.ascending:
            reactioni, reactionj = reactionj, reactioni
            self.momenti, self.momentj = -self.momentj, -self.momenti
        self.reactioni = globalize_perpendicular_load(reactioni, element.angle)
        self.reactionj = globalize_perpendicular_load(reactionj, element.angle)


def global_matrixFrame(elements:List[ElementFrame]):
    gmatrix = np.zeros((Node3DoF.code*3, Node3DoF.code*3))
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


def load_elements(elements:List[ElementFrame], loads:List[List[PerpendicularLoad]]):
    assert len(loads) == len(elements), "Loads must equal elements"
    
    load_matrix = np.zeros((Node3DoF.code*3, 1))
    for i in range(len(loads)):
        if loads[i]:
            for load in loads[i]:
                load.apply(elements[i])
                inode = elements[i].get_nodes()[0].get_code()
                jnode = elements[i].get_nodes()[1].get_code()
                load_matrix[inode*3] += load.reactioni[0]
                load_matrix[inode*3+1] += load.reactioni[1]
                load_matrix[inode*3+2] += load.momenti
                load_matrix[jnode*3] += load.reactionj[0]
                load_matrix[jnode*3+1] += load.reactionj[1]
                load_matrix[jnode*3+2] += load.momentj
    return load_matrix


if __name__ == "__main__":
    #Example of usage:
    #Start by defining the nodes (x, y, name and if there is any limitation to their movement [X, Y| X, Y, Z | X,Y, rotation])
    #This depends on the elements being modeled
    
    nodes = [Node3DoF(0, 0, "A", (True, True, False)), Node3DoF(480, 0, "B"), Node3DoF(480, -160, "C"), Node3DoF(0, -400, "D", (True, False, False))]
    
    #Definition of Frame Elements (2D), young modulus, area, second moment of area and the nodes that make up each element. Node i, then j.
    elements = [ElementFrame(200, 480, 77e3, (nodes[0], nodes[1])), ElementFrame(200, 480, 77e3, (nodes[1], nodes[2])), ElementFrame(200, 480, 77e3, (nodes[2], nodes[3]))]
    
    #Create global stiffness matrix from the elements just created
    matrizg = global_matrixFrame(elements)

    #Create load matrix
    #Element 1 has a uniform load of 0.2N/mm and a triangular load of 0.2N/mm with a negative slope so it is especified as no ascendig
    #Element 2 has no load and has a point load 386.66mm from the node i of 250N
    load_matrix = load_elements(elements, [(UniformLoad(0.2), TriangleLoad(0.2,ascending=False)), (), (PointLoad(250, 386.66),)])

    #Obtain displacement and reaction, interpreting 3D dimention from the nodes as rotation freedom/restriction instead of movement in Z direction for this particular problem
    dis, reac = solve_matrix3D(matrizg, load_matrix, nodes)
    
    
    #Printing results
    coord = ["x", "y", "rot."]
    units = ["mm", "mm", "rad"]
    units2 = ["N", "N", "N-mm"]
    print("-"*22+"Displacements"+"-"*22)
    for i in range(1, len(dis)+1):
        print(f"Displacement from the node {nodes[(i-1)//3].name} en {coord[(i-1)%3]} es:\t\t\t{dis[i-1][0]} {units[(i-1)%3]}")
    print("-"*60)
    print("-"*25+"Reacciones"+"-"*25)
    for i in range(1, len(reac)+1):
        print(f"Reaction form node {nodes[(i-1)//3].name} en {coord[(i-1)%3]} es:\t\t\t{reac[i-1][0]} {units2[(i-1)%3]}")
    
    
    
    
    
    
