from examples.rings import Z
from algorithms.gcd_dfu import gcd_dfu


#En primer lugar, cogemos el UFD que será la base de nuestro anillo de polinomios
ZX = Z["X"]

#TODO: probar con polinomios que no tengan raíces comunes



#Polinomios de ejemplo
res = ZX.build([-1,0,1,1])
f = ZX.build([-5,-2,9,4,-3,-2,4,2,1])
g = ZX.build([6,7,-11,-13,-2,5])


#Imprimimos el GCD:
print("Test 1:")
print("f: ",f)
print("g: ",g)
print(gcd_dfu(f, g))
print(gcd_dfu(-f, g))
print(gcd_dfu(f, -g))
print(gcd_dfu(-f, -g))
if gcd_dfu(f, g) == res:
    print("Correcto\n")
else:
    print("Incorrecto\n")


#(x^2+8x+3)(x-17) = x^3-9x^2-133x-51
#(x^2+8x+3)(x-17)(x^3-2) = x^6-9x^5-133x^4-53x^3+18x^2+266x+102
res = ZX.build([-51,-133,-9,1])
a = ZX.build([-51,-133,-9,1])
b = ZX.build([102,266,18,-53,-133,-9,1])

print("Test 2:")
print("f: ",a)
print("g: ",b)
gcd = gcd_dfu(a, b)
print("GCD: ", gcd)
if gcd == res:
    print("Correcto\n")
else:
    print("Incorrecto\n")


a = a * ZX.build([7,1])
print("Test 3:")
print("f: ",a)
print("g: ",b)
gcd = gcd_dfu(a, b)
print("GCD: ", gcd)
if gcd == res:
    print("Correcto\n")
else:
    print("Incorrecto\n")
