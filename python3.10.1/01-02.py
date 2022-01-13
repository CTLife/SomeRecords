print('##1.start####################################################################################\n')
print("Hello World.")
print("Hello\n\n Wo\trld.")
print('Hello\n\n Wo\t\t\trld.')

a = 5  # 定义一个数值类型变量
print(a)
print(type(a))
a = 'Hello, Charlie'  # 重新将字符串值赋值给a变量
print(a)
print(type(a))

ac1 = 3 + 0.2j
print(ac1)
print(type(ac1)) # 输出 complex类型
ac2 = 4 - 0.1j
print(ac2)
print(ac1 + ac2) # # 复数运行, 输出 (7+0.1j)

import cmath  # 导入cmatch模块
ac3 = cmath.sqrt(-1)  # sqrt()是cmath模块下的函数，用于计算平方根
print(ac3) # 输出 1j
print(1j**2)
print((-1)**0.5)
print('##1.end####################################################################################\n\n\n\n\n')





print('##2.start####################################################################################\n')
# 这是一行简单的注释
print("Hello World!")
'''
这里面的内容全部是多行注释.
Python语言真的很简单.
'''
# print("这行代码被注释了，将不会被编译、执行!")
"""
这是用三个双引号括起来的多行注释。
Python同样是允许的。
"""

af1 = 5.2345556
print("af1的值为:", af1) # 输出af1的值

af2 = 25.2345
print("af2的类型为:", type(af2))

f1 = 5.12e2
print("f1的值为:", f1)
f2 = 5e3
print("f2的值为:", f2)
print("f2的类型为:", type(f2)) # 看到类型为float


# 以0x或0X开头的整数数值是十六进制的整数
hex_value1 = 0x13
hex_value2 = 0XaF
print("hexValue1的值为：", hex_value1)
print("hexValue2的值为：", hex_value2)

# 以0b或0B开头的整数数值是二进制的整数
bin_val = 0b111
print('bin_val的值为：', bin_val)
bin_val = 0B101
print('bin_val的值为：', bin_val)

# 以0o或0O开头的整数数值是8进制的整数
oct_val = 0o54
print('oct_val的值为：', oct_val)
oct_val = 0O17
print('oct_val的值为：', oct_val)

# 在数值中使用下画线
one_million = 1_000_000
print(one_million)
price = 234_234_234 # price实际的值为234234234
android = 1234_1234 # android实际的值为12341234
print(price+android)
print(2_1+1_5)
print('##2.end####################################################################################\n\n\n\n\n')





print('##3.start####################################################################################\n')
user_name = 'Charlie'
user_age = 8
print("读者名:" , user_name, "年龄:", user_age)            # 同时输出多个变量和字符串
print("读者名:" , user_name, "年龄:", user_age, sep='|')   # 同时输出多个变量和字符串，指定分隔符
print("读者名:" , user_name, "年龄:", user_age, sep='\n')  # 同时输出多个变量和字符串，指定分隔符
print("读者名:" , user_name, "年龄:", user_age, sep='##')  # 同时输出多个变量和字符串，指定分隔符

# 指定end参数，指定输出之后不再换行
print(40, '\t', end="")
print(50, '\t', end="")
print(60, '\t', end="")
print(70, '\t', end="\n\n\n")


f = open("poem.txt", "w") # 打开文件以便写入
print('沧海月明珠有泪', file=f)
print('蓝田日暖玉生烟', file=f)
f.close()

import keyword
print( keyword.kwlist )

a = 56  # 定义变量a，赋值为56
print(a)
a = 9999999999999999999999  # 为a赋值一个大整数
print(a)
print(type(a))  # type()函数用于返回变量的类型
a = None
print(a)
print(type(a))
print('##3.end####################################################################################\n\n\n\n\n')





print('##4.start####################################################################################\n')
b1 = bytes()   # 创建一个空的bytes
b2 = b''       # 创建一个空的bytes值
b3 = b'hello'  # 通过b前缀指定hello是bytes类型的值
print(b1)
print(b2)
print(b3)
print(b3[0])
print(b3[2:4])

b4 = bytes('我爱Python编程',encoding='utf-8')  # 调用bytes方法将字符串转成bytes对象
print(b4)
b5 = "学习Python很有趣".encode('utf-8')  # 利用字符串的encode()方法编码成bytes，默认使用utf-8字符集
print(b5)
st = b5.decode('utf-8')  # 将bytes对象解码成字符串，默认使用utf-8进行解码。
print(st) # 学习Python很有趣

msg = input("请输入你的值：")
print(type(msg))
print(msg)
print( int(msg)+5 )  ## print( msg+5 ) is wrong
print( float(msg)+5 )

s = '''"Let's go fishing", said Mary.
"OK, Let's go", said her brother.
they walked to a lake.
'''
print(s)

s2 = 'The quick brown fox \
jumps over the lazy dog.'
print(s2)

num = 20 + 3 / 4 + \
    2 * 3
print(num)

s1 = r'G:\publish\codes\02\2.4'
print(s1)
s2 = r'"Let\'s go", said Charlie.'  # 原始字符串包含的引号，同样需要转义
print(s2)
s3 = r'Good Morning' '\\'
print(s3)
print('##4.end####################################################################################\n\n\n\n\n')





print('##5.start####################################################################################\n')
s1 = "Hello,"  'Charlie'
print(s1)

s2 = "Python "
s3 = "is Funny"
# 使用+拼接字符串
s4 = s2 + s3
print(s4) 

s5 = '123' '456' "789"
print(s5)

str1 = 'Charlie'
str2 = "疯狂软件教育"
print(str1)
print(str2)
str3 = "I'm a coder"  #str3 = 'I'm a coder'  is wrong.
str4 = '"Spring is here, let us jam!", said woodchuck.'
str5 = '"we are scared, Let\'s hide in the shade", says the bird'
print(str5)

s1 = "这本书的价格是："
p = 99.8
#print(s1 + p) # 字符串直接拼接数值，程序报错
print(s1 + str(p))  # 使用str()将数值转换成字符串
print(s1 + repr(p)) # 使用repr()将数值转换成字符串

st = "I will play my fife."
print(st)
print(repr(st)) 
print('##5.end####################################################################################\n\n\n\n\n')





print('##6.start####################################################################################\n')
a = 'our domain is crazyit.org'
print(a.title())  # 每个单词首字母大写
print(a.lower())  # All字母小写
print(a.upper())  # All字母大写

s = 'crazyit.org is very good'
print(s[2])       # 获取s中索引2处的字符 # 输出a
print(s[-4])      # 获取s中从右边开始，索引4处的字符 # 输出g
print(s[3: 5])    # 获取s中从索引3处到索引5处（不包含）的子串  # 输出zy
print(s[3: -5])   # 获取s中从索引3处到倒数第5个字符的子串  # 输出zyit.org is very
print(s[-6: -3])  # 获取s中从倒数第6个字符到倒数第3个字符的子串  # 输出y g
print(s[5: ])     # 获取s中从索引5处到结束的子串  # 输出it.org is very good
print(s[-6: ])    # 获取s中从倒数第6个字符到结束的子串  # 输出y good
print(s[: 5])     # 获取s中从开始到索引5处的子串  # 输出crazy
print(s[: -6])    # 获取s中从开始到倒数第6个字符的子串  #输出crazyit.org is ver
print('very' in s) # 判断s是否包含'very'子串  # True
print('fkit' in s) # False
print(len(s))      # 输出s的长度# 24
print(len('test')) # 输出'test'的长度  # 4
print(max(s))      # 输出s字符串中最大的字符  # z
print(min(s))      # 输出s字符串中最大的字符   # 空格

s = 'Hello\nCharlie\nGood\nMorning.'
print(s)
s2 = '商品名\t\t单价\t\t数量\t\t总价'
s3 = '疯狂Java讲义\t108\t\t2\t\t316'
print(s2)
print(s3)
print('##6.end####################################################################################\n\n\n\n\n')





print('##7.start####################################################################################\n')
price = 108
print("the book's price is %x" % price)
print("the book's price is %f" % price)
print("the book's price is %d" % price)
print("the book's price is %i" % price)
print("#######################\n\n")

user = "Charli"
age = 8
print("%s is a %s years old boy" % (user , age))  # 格式化字符串有两个占位符，第三部分提供2个变量
num = -28
print("num is: %6i" % num)
print("num is: %6d" % num)
print("num is: %6o" % num)
print("num is: %6x" % num)
print("num is: %6X" % num)
print("num is: %6s" % num)
print("#######################\n\n")

num2 = 30
print("num2 is: %06d" % num2)  # 最小宽度为0，左边补0
print("num2 is: %+06d" % num2) # 最小宽度为6，左边补0，总带上符号
print("num2 is: %-6d" % num2)  # 最小宽度为6，右对齐

my_value = 3.001415926535
print("my_value is: %8.3f" % my_value)    # 最小宽度为8，小数点后保留3位
print("my_value is: %08.3f" % my_value)   # 最小宽度为8，小数点后保留3位，左边补0
print("my_value is: %+08.3f" % my_value)  # 最小宽度为8，小数点后保留3位，左边补0，始终带符号

my_value = -3.001415926535
print("my_value is: %8.3f" % my_value)    # 最小宽度为8，小数点后保留3位
print("my_value is: %08.3f" % my_value)   # 最小宽度为8，小数点后保留3位，左边补0
print("my_value is: %+08.3f" % my_value)  # 最小宽度为8，小数点后保留3位，左边补0，始终带符号
print("#######################\n\n")

the_name = "Charlie"
print("the name is: %.3s" % the_name)   # 只保留3个字符  # 输出Cha
print("the name is: %10.2s" % the_name) # 只保留2个字符，最小宽度10
print('##7.end####################################################################################\n\n\n\n\n')





print('##8.start####################################################################################\n')
s = 'crazyit.org is a good site'
print(s.startswith('crazyit'))  # 判断s是否以crazyit开头
print(s.endswith('site'))       # 判断s是否以site结尾
print(s.find('org'))            # 查找s中'org'的出现位置  # 8
print(s.index('org'))           # 查找s中'org'的出现位置  # 8
print(s.find('org', 9))         # 从索引为9处开始查找'org'的出现位置 # -1
#print(s.index('org', 9))       # 从索引为9处开始查找'org'的出现位置  # 引发错误

print(s.replace('it', 'xxxx'))        # 将字符串中所有it替换成xxxx
print(s.replace('it', 'xxxx', 1))     # 将字符串中1个it替换成xxxx
table = {97: 945, 98: 946, 116: 964}  # 定义替换表：97（a）->945（α）,98（b）->945（β）,116（t）->964（τ）,
print(s.translate(table))             # crαzyit.org is α good site
table1 = str.maketrans('is', '89')
print(s.translate(table1))

s = 'crazyit.org is a good site'
print(s.split())         # 使用空白对字符串进行分割  # 输出 ['crazyit.org', 'is', 'a', 'good', 'site']
print(s.split(None, 2))  # 使用空白对字符串进行分割,最多只分割前2个单词  # 输出 ['crazyit.org', 'is', 'a good site']
print(s.split('.'))      # 使用点进行分割 # 输出 ['crazyit', 'org is a good site']
mylist = s.split()
print('/'.join(mylist))  # 使用'/'为分割符，将mylist连接成字符串  # 输出 crazyit.org/is/a/good/site
print(','.join(mylist))  # 使用','为分割符，将mylist连接成字符串  # 输出 crazyit.org,is,a,good,site

s = '  this is a puppy            '
print(s.lstrip())   # 删除左边的空白
print(s.rstrip())   # 删除右边的空白
print(s.strip())    # 删除两边的空白
print(s)            # 再次输出s，将会看到s并没有改变

s2 = 'i think it is a scarecrow'
print(s2.lstrip('itow'))  # 删除左边的i、t、o、w字符
print(s2.rstrip('itow'))  # 删除右边的i、t、o、w字符
print(s2.strip('itow'))   # 删除两边的i、t、o、w字符
print('##8.end####################################################################################\n\n\n\n\n')





print('##9.start####################################################################################\n')
a = 5.2
b = 3.1
the_sum = a + b
print("the_sum的值为：", the_sum)  # sum的值为8.3

s1 = 'Hello, '
s2 = 'Charlie'
print(s1 + s2)  # 使用+连接两个字符串
print('1' + '2' + '3')

c = 5.2
d = 3.1
sub = c - d
print("sub的值为：", sub)  # sub的值为2.1

e = 5.2
f = 3.1
multiply = e * f
print("multiply的值为：", multiply)  # multiply的值为16.12

s3 = 'crazyit '
print(s3 * 5)  # 使用*将5个字符串连接起来

print("19/4的结果是:", 19/4)
print("19//4的结果是:", 19//4)
aa = 5.2
bb = 3.1
print("aa/bb的值是:", aa / bb)   # aa / bb的值将是1.67741935483871
print("aa//bb的值是:", aa // bb) # aa // bb值将是1.0

print('5的2次方：', 5 ** 2)            # 25
print('4的3次方：', 4 ** 3)           # 64
print('4的开平方：', 4 ** 0.5)        # 2.0
print('27的开3次方：', 27 ** (1 / 3)) # 3.0

st = "Python"    # 为变量st赋值为Python
pi = 3.14        # 为变量pi赋值为3.14
visited  = True  # 为变量visited赋值为True
print(st)
print(pi)
print(visited)

st2 = st  # 将变量st的值赋给st2
print(st2)

a = b = c = 20
print(a)
print(b)
print(c)

d1 = 12.34
d2 = d1 + 5  # 将表达式的值赋给d2
# 输出d2的值
print("d2的值为：%g" % d2 ) # 17.34
print("d2的值为：%f" % d2 ) # 17.34
print("d2的值为：%.3f" % d2 ) # 17.34
print('##9.end####################################################################################\n\n\n\n\n')








print('##10.start####################################################################################\n')
print(5 & 9)  # 将输出1
print(5 | 9)  # 将输出13

a = -5
print( ~a)    # 将输出4
print(5 ^ 9)  # 将输出12

print(5 << 2)  # 输出20
print(-5 << 2) # 输出-20

b = -5
print(b >> 2)  # 输出-2

bookName = "疯狂Python"
price = 79
version = "正式版"
if bookName.endswith('Python') and (price < 50 or version == "正式版") :
    print("打算购买这本Python图书")
else:
    print("不购买！")

print("5是否大于 4：", 5 > 4)                         # 输出True
print("3的4次方是否大于等于90.0：", 3 ** 4 >= 90)      # 输出False
print("20是否大于等于20.0：", 20 >= 20.0)             # 输出True
print("5和5.0是否相等：", 5 == 5.0)                  # 输出True
print("True和False是否相等：", True == False)        # 输出False


print("1和True是否相等：", 1 == True)     # 输出True
print("0和False是否相等：", 0 == False)   # 输出True
print(True + False)  # 输出1
print(False - True)  # 输出-1
print('##10.end####################################################################################\n\n\n\n\n')








print('##11.start####################################################################################\n')
import time
a = time.gmtime()  # 获取当前时间
b = time.gmtime()
print(a == b) # a和b两个时间相等，输出True
print(a is b) # a和b不是同一个对象，输出False
print(id(a))
print(id(b))

a = 'abcdefghijklmn'
print(a[2:8:3])   # 获取索引2到索引8的子串，步长为3  # 输出cf
print(a[2:8:2])   # 获取索引2到索引8的子串，步长为2  # 输出ceg
print(a[2:8:1])
 
s = 'crazyit.org'
print('it' in s)       # True
print('it' not in s)   # False
print('fkit' in s)     # False
print('fkit' not in s) # True

print(not False)            # 直接对False求非运算，将返回True
print(5 > 3 and 20.0 > 10)  # 5>3返回True，20.0大于10，因此结果返回True
print(4 >= 5 or "c" > "a")  # 4>=5返回False，"c">"a"返回True。求或后返回True
print('##11.end####################################################################################\n\n\n\n\n')










print('##12.start####################################################################################\n')
print("5%3的值为：", 5 % 3)             # 输出2
print("5.2%3.1的值为：", 5.2 % 3.1)     # 输出2.1
print("-5.2%-3.1的值为：", -5.2 % -3.1) # 输出-2.1
print("5.2%-2.9的值为：", 5.2 % -2.9)   # 输出-0.6
print("5.2%-1.5的值为：", 5.2 % -1.5)   # 输出-0.8
print("-5.2%1.5的值为：", -5.2 % 1.5)   # 输出0.8
#print("5对0.0求余的结果是:", 5 % 0.0)   # 导致错误

print("5/3的值为：", 5 // 3)              
print("5.2//3.1的值为：", 5.2 // 3.1)      
print("-5.2//-3.1的值为：", -5.2 // -3.1)  
print("5.2//-2.9的值为：", 5.2 // -2.9)    
print("5.2//-1.5的值为：", 5.2 // -1.5)   
print("-5.2//1.5的值为：", -5.2 // 1.5)    
 
x = -5.0  # 定义变量x，其值为-5.0
x = -x    # 将x求负，其值变成5.0
print(x)

y = -5.0  # 定义变量y，其值为-5.0
y = +y    # y值依然是-5.0
print(y)

a = 5
b = 3
st = "a大于b" if a > b else  "a不大于b" 
print(st)  # 输出"a大于b"

print("a大于b") if a > b else print("a不大于b") # 输出"a大于b"

# 第一个返回值部分使用两条语句，逗号隔开
st = print("crazyit"), 'a大于b' if a > b else  "a不大于b" 
print(st)

# 第一个返回值部分使用两条语句，分号隔开
st = print("crazyit"); x = 20 if a > b else  "a不大于b" 
print(st)
print(x)

c = 5
d = 5
# 下面将输出c等于d
print("c大于d") if c > d else (print("c小于d") if c < d else print("c等于d"))
print('##12.end####################################################################################\n\n\n\n\n')

















