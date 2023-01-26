import random

open("test_A.fa", "x")
f = open("test_A.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    f.write("AAAAAAAA\n")

open("test_G.fa", "x")
f = open("test_G.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    f.write("GGGGGGGG\n")

open("test_T.fa", "x")
f = open("test_T.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    f.write("TTTTTTTT\n")

open("test_C.fa", "x")
f = open("test_C.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    f.write("CCCCCCCC\n")

open("test_AGTC.fa", "x")
f = open("test_ATGC.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    f.write("AGTCAGTC\n")

open("test_rand_1.fa", "x")
f = open("test_rand_1.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n)+ "-" + str(n+8) + "(-)\n" )
    n = n + 8
    for k in range(8):
        f.write(random.choice("ATGC"))
    f.write("\n")

open("test_rand_2.fa", "x")
f = open("test_rand_2.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    for k in range(8):
        f.write(random.choice("ATGC"))
    f.write("\n")

open("test_rand_3.fa", "x")
f = open("test_rand_3.fa", "a")
n = 0
for i in range(1000):
    f.write(">Chr1:" + str(n) + "-" + str(n+8) + "(-)\n" )
    n = n + 8
    for k in range(8):
        f.write(random.choice("ATGC"))
    f.write("\n")