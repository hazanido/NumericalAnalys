import time

input = int(input('Enter id:'))

random = (((time.localtime().tm_hour+time.localtime().tm_min)*time.localtime().tm_sec)//13) + input

print((random % 12) + 19)



