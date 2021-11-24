#week 1 lab session 

x = 5
y = 11
z = 18
xyz = x * y * z
print(xyz)
ln_xyz = log(xyz)
print(ln_xyz)
check = ln_xyz > x
print(check)

fr_ln = seq(from= ln_xyz, to= xyz,by= 10)
print(fr_ln)


length(fr_ln)

new_var = c(fr_ln,x,y,z)
print(length(new_var))

