# Gates commonly used for yeast and E. coli live cell cytometry
# TODO: test each gate, make sure it's the same as in R 

from fcm import PolyGate

# yeastGate: gates on SSC-A vs. FSC-A for yeast-sized objects on Accuri C6
yeastGate = PolyGate([[4e5,1e4],
                      [1e7,1e4],
                      [1e7,2.3e6],
                      [4e5,6e4]],
                     [0,1],
                     'yeastGate')

dipsingletGate = PolyGate([[7.5e5,9e5],
                           [1.3e6,1.6e6],
                           [1.8e6,2.6e6],
                           [1.5e6,3e6],
                           [6e5,1.5e6]],
                          [6,0],
                          'dipsingletGate')

dipdoubletGate = PolyGate([[1e6,8e5],
                           [1.7e6,1.25e6],
                           [2.3e6,1.7e6],
                           [2.2e5,2e6],
                           [2e6,2.2e6],
                           [8e5,8.5e5]],
                          [6,0],
                          'dipdoubletGate')

hapsingletGate = PolyGate([[5e5,8e5],
                           [8e5,1.05e6],
                           [1.15e6,1.5e6],
                           [1e6,1.8e6],
                           [5e5,1e6]],
                          [6,0],
                          'hapsingletGate')

hapdoubletGate = PolyGate([[6.5e5,5.75e5],
                           [1.15e6,9e5],
                           [1.5e6,1.3e6],
                           [1.4e6,1.4e6],
                           [1.2e6,1.5e6],
                           [5e5,6.5e5]],
                          [6,0],
                          'hapdoubletGate')



