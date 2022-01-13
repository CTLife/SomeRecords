print('##1.start####################################################################################\n')
my_list = ['crazyit', 20, 'Python']  # 使用方括号定义列表
print(my_list)
my_tuple = ('crazyit', 20, 'Python') # 使用圆括号定义元组
print(my_tuple)

my_list[1]=88
#my_tuple[1]=88
print(my_list)
print(my_tuple)

a_tuple = (20, 10, -2, 15.2, 102, 50)  # 元素都是数值的元组
print(max(a_tuple))   # 计算最大值 # 102
print(min(a_tuple))   # 计算最小值 # -2
print(len(a_tuple))   # 计算长度  # 6

b_list = ['crazyit', 'fkit', 'Python', 'Kotlin']  # 元素都是字符串的列表
print(max(b_list))   # 计算最大值（依次比较每个字符的ASCII码值，先比较第一个字符，若相同，继续比较第二个字符，以此类推） # fkit（26个小写字母的ASCII码为97~122）
print(min(b_list))   # 计算最小值 # Kotlin （26个大写字母的ASCII码为65~90）
print(len(b_list))   # 计算长度 # 4

a_tuple = ('crazyit', 20, 5.6, 'fkit', -17)
print(a_tuple)
print(a_tuple[0])   # 访问第1个元素  # crazyit
print(a_tuple[1])   # 访问第2个元素  # 20
print(a_tuple[-1])  # 访问倒数第1个元素 # -17
print(a_tuple[-2])  # 访问倒数第2个元素 # -fkit

a_list = ['crazyit', 20, 5.6, 'fkit', -17]
print(a_list)
print(a_list[0])   # 访问第1个元素  # crazyit
print(a_list[1])   # 访问第2个元素  # 20
print(a_list[-1])  # 访问倒数第1个元素 # -17
print(a_list[-2])  # 访问倒数第2个元素 # -fkit

a_tuple = ('crazyit' , 20, -1.2)
print(20 in a_tuple)          # True
print(1.2 in a_tuple)         # False
print('fkit' not in a_tuple)  # True
print('##1.end####################################################################################\n\n\n\n\n')




print('##2.start####################################################################################\n')
a_tuple = ('crazyit' , 20)
mul_tuple = a_tuple * 3  # 执行乘法
print(mul_tuple)         # ('crazyit', 20, 'crazyit', 20, 'crazyit', 20)

a_list = [30, 'Python', 2]
mul_list = a_list * 3
print(mul_list) # [30, 'Python', 2, 30, 'Python', 2, 30, 'Python', 2]

a_tuple = ('crazyit' , 20, -1.2)
b_tuple = (127, 'crazyit', 'fkit', 3.33)
sum_tuple = a_tuple + b_tuple  # 计算元组相加
print(sum_tuple) # ('crazyit', 20, -1.2, 127, 'crazyit', 'fkit', 3.33)
print(a_tuple)  # a_tuple并没有改变
print(b_tuple)  # b_tuple并没有改变
print(a_tuple + (-20 , -30)) # 两个元组相加 # ('crazyit', 20, -1.2, -20, -30)

# 下面代码报错：元组和列表不能直接相加
#print(a_tuple + [-20 , -30])
a_list = [20, 30, 50, 100]
b_list = ['a', 'b', 'c']
sum_list = a_list + b_list # 计算列表相加
print(sum_list)           # [20, 30, 50, 100, 'a', 'b', 'c']
print(a_list + ['fkit']) # [20, 30, 50, 100, 'fkit']
print(sum_list[2]+50)

# 序列封包：将10、20、30封装成元组后赋值给vals
vals = 10, 20, 30
print(vals)        # (10, 20, 30)
print(type(vals)) # <class 'tuple'>
print(vals[1])   # 20

a_tuple = tuple(range(1, 10, 2))
# 序列解包: 将a_tuple元组的各元素依次赋值给a、b、c、d、e变量
a, b, c, d, e = a_tuple
print(a, b, c, d, e) # 1 3 5 7 9

a_list = ['fkit', 'crazyit']
# 序列解包: 将a_list序列的各元素依次赋值给a_str、b_str变量
a_str, b_str = a_list
print(a_str, b_str) # fkit crazyit


x, y, z = 10, 20, 30  # 将10、20、30依次赋值给x、y、z
print(x, y, z)       # 10 20 30
x, y, z = y, z, x   # 将y,z, x依次赋值给x、y、z
print(x, y, z)     # 20 30 10

first, second, *rest = range(10)  # first、second保存前2个元素，rest列表包含剩下的元素
print(first)  # 0
print(second) # 1
print(rest)   # [2, 3, 4, 5, 6, 7, 8, 9]

*begin, last = range(10)  # last保存最后一个元素，begin保存前面剩下的元素
print(begin) # [0, 1, 2, 3, 4, 5, 6, 7, 8]
print(last) # 9

first, *middle, last = range(10)  # first保存第一个元素，last保存最后一个元素，middle保存中间剩下的元素
print(first)   # 0
print(middle) # [0, 1, 2, 3, 4, 5, 6, 7, 8]
print(last)  # 9
print('##2.end####################################################################################\n\n\n\n\n')




print('##3.start####################################################################################\n')

a_tuple = ('crazyit', 20, 5.6, 'fkit', -17)
print(a_tuple[1: 3])      # 访问从第2个到倒数第4个（不包含）所有元素  # (20, 5.6)
print(a_tuple[-3: -1])    # 访问从倒数第3个到倒数第1个（不包含）所有元素  # (5.6, 'fkit')
print(a_tuple[1: -2])     # 访问从第2个到倒数第2个（不包含）所有元素  # (20, 5.6)
print(a_tuple[-3: 4])     # 访问从倒数第3个到第5个（不包含）所有元素  # (5.6, 'fkit')

b_tuple = (1, 2, 3, 4, 5, 6, 7, 8, 9)
print(b_tuple[2: 8: 2])   # 访问从第3个到第9个（不包含）、间隔为2的所有元素  # (3, 5, 7)
print(b_tuple[2: 8: 3])   # 访问从第3个到第9个（不包含）、间隔为3的所有元素  # (3, 6)
print(b_tuple[2: -2: 2])  # 访问从第3个到倒数第2个（不包含）、间隔为3的所有元素 # (3, 5, 7)
print(b_tuple[2: 8: 1])

a_list = ['crazyit', 20, 5.6, 'fkit', -17]
print(a_list[1: 3])    # 访问从第2个到倒数第4个（不包含）所有元素  # (20, 5.6)
print(a_list[-3: -1])  # 访问从倒数第3个到倒数第1个（不包含）所有元素  # (5.6, 'fkit')
print(a_list[1: -2])   # 访问从第2个到倒数第2个（不包含）所有元素 # (20, 5.6)
print(a_list[-3: 4])   # 访问从倒数第3个到第5个（不包含）所有元素 # (5.6, 'fkit')

b_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
print(b_list[2: 8: 2])  # 访问从第3个到第9个（不包含）、间隔为2的所有元素 # (3, 5, 7)
print(b_list[2: 8: 3])  # 访问从第3个到第9个（不包含）、间隔为3的所有元素 # (3, 6)
print(b_list[2: -2: 2]) # 访问从第3个到倒数第2个（不包含）、间隔为3的所有元素 # (3, 5, 7)
print(b_list[2: 8: 1])


# 同时对元组使用加法、乘法
order_endings = ('st', 'nd', 'rd')\
    + ('th',) * 17 + ('st', 'nd', 'rd')\
    + ('th',) * 7 + ('st',)

print(order_endings) # 将会看到st、nd、rd、17个th、st、nd、rd、7个th、st
day = input("输入日期(1-31)：")
day_int = int(day)  # 将字符串转成整数
print(day + order_endings[day_int - 1])


a1=(1)
a2=(1,)
b1=[1]
b2=[1,]

print(a1, type(a1))
print(a2, type(a2))
print(b1, type(b1))
print(b2, type(b2))
print('##3.end####################################################################################\n\n\n\n\n')




print('##4.start####################################################################################\n')
a_list = ['crazyit', 20, -2]
a_list.append('fkit')  # 追加元素
print(a_list)          # ['crazyit', 20, -2, 'fkit']

a_tuple = (3.4, 5.6)
a_list.append(a_tuple) # 追加元组，元组被当成一个元素
print(a_list)          # ['crazyit', 20, -2, 'fkit', (3.4, 5.6)]
print(a_list[3])
print(a_list[4])


a_list.append(['a', 'b'])  # 追加列表，列表被当成一个元素
print(a_list)              # ['crazyit', 20, -2, 'fkit', (3.4, 5.6), ['a', 'b']]
print(a_list[3])
print(a_list[5])


a_list = ['crazyit', 20, -2]
a_list[-1:-1] = 'fkit' # 追加元素
print(a_list) 

a_tuple = (3.4, 5.6)
a_list[-1:-1] = a_tuple
print(a_list)  
a_list[0:0] = a_tuple
print(a_list)  
 

# 追加列表，列表被当成一个元素
a_list.append(['a', 'b'])
print(a_list) # ['crazyit', 20, -2, 'fkit', (3.4, 5.6), ['a', 'b']]
print(a_list[3])
print(a_list[5])




b_list = ['a', 30]
b_list.extend((-2, 3.1))  ## 追加元组中的所有元素
print(b_list)             # ['a', 30, -2, 3.1]

b_list.extend(['C', 'R', 'A'])  # 追加列表中的所有元素
print(b_list)                   # ['a', 30, -2, 3.1, 'C', 'R', 'A']

b_list.extend(range(97, 100))  # 追加区间中的所有元素
print(b_list)                  # ['a', 30, -2, 3.1, 'C', 'R', 'A', 97, 98, 99]

c_list = list(range(1, 6))
print(c_list)                 # [1, 2, 3, 4, 5]
c_list.insert(3, 'CRAZY' )    # 在索引3处插入字符串
print(c_list)                 # [1, 2, 3, 'CRAZY', 4, 5]
c_list.insert(3, tuple('crazy'))  # 在索引3处插入元组，元组被当成一个元素
print(c_list)                     # [1, 2, 3, ('c', 'r', 'a', 'z', 'y'), 'CRAZY', 4, 5]

# assign lists
list_1 = [1, 2, 3]
list_2 = [1, 2, 3]
list_3 = [1, 2, 3]
 
a = [2, 3]
# use methods
list_1.append(a)
list_2.insert(3, a)
list_3.extend(a)
 
# display lists
print(list_1)
print(list_2)
print(list_3)

# Python3 code to demonstrate to insert one list in another using list slicing.
  
# initializing lists 
test_list = [4, 5, 6, 3, 9]
insert_list = [2, 3]
  
# initializing position
pos = 2
  
# printing original list
print ("The original list is : " + str(test_list))
  
# printing insert list 
print ("The list to be inserted is : " + str(insert_list))
  
# using list slicing
# to insert one list in another
test_list[pos:pos] = insert_list
  
# printing result 
print ("The list after insertion is : " +  str(test_list))


a_list = [2, 30, 'a', [5, 30], 30]
print(a_list.count(30))         # 计算列表中30的出现次数  # 2
print(a_list.count([5, 30]))    # 计算列表中[5, 30]的出现次数  # 1

dir(list)

a = range(1,11)
print(type(a))
print( a )
print( list(a) )
print('##4.end####################################################################################\n\n\n\n\n')




print('##5.start####################################################################################\n')
a_list = ['crazyit', 20, -2.4, (3, 4), 'fkit']
del a_list[2]      # 删除第3个元素
print(a_list)      # ['crazyit', 20, (3, 4), 'fkit']
del a_list[1: 3]   # 删除第2个到第4个（不包含）元素
print(a_list)      # ['crazyit', 'fkit']

b_list = list(range(1, 10))
del b_list[2: -2: 2]  # 删除第3个到倒数第2个（不包含）元素，间隔为2
print(b_list)         # [1, 2, 4, 6, 8, 9]
del b_list[2: 4]      # 删除第3个到第5个（不包含）元素
print(b_list)         # [1, 2, 8, 9]

name = 'crazyit'
print(name)  # crazyit
del name     # 删除name变量
#print(name) # NameError

c_list = [20, 'crazyit', 30, -4, 'crazyit', 3.4]
c_list.remove(30)          # 删除第一次找到的30
print(c_list)              # [20, 'crazyit', -4, 'crazyit', 3.4]
c_list.remove('crazyit')   # 删除第一次找到的'crazyit'
print(c_list)              # [20, -4, 'crazyit', 3.4]

c_list.clear()
print(c_list) # []

a_list = [2, 30, 'a', 'b', 'crazyit', 30]
print(a_list.index(30))       # 定位元素30的出现位置  # 1
print(a_list.index(30, 2))    # 从索引2处开始、定位元素30的出现位置  # 5
#print(a_list.index(30, 2, 4)) # 从索引2处到索引4处之间定位元素30的出现位置，找不到该元素 # ValueError

a_tuple = ('crazyit', 20, -1.2)
print(type(a_tuple))
a_list = list(a_tuple)  # 将元组转换成列表
print(a_list)
print(type(a_list))


a_range = range(1, 5)   # 使用range()函数创建区间（range）对象
print(a_range)          # range(1, 5)
b_list = list(a_range)  # 将区间转换成列表
print(b_list)           #[1, 2, 3, 4]

c_list = list(range(4, 20, 3))   # 创建区间时还指定步长
print(c_list)                    # [4, 7, 10, 13, 16, 19]

a_list = list(range(1, 8))
a_list.reverse()   # 将a_list列表元素反转
print(a_list)      # [7, 6, 5, 4, 3, 2, 1]
print('##5.end####################################################################################\n\n\n\n\n')






print('##6.start####################################################################################\n')
a_list = [3, 4, -2, -30, 14, 9.3, 3.4]
a_list.sort()  # 对列表元素排序
print(a_list)  #[-30, -2, 3, 3.4, 4, 9.3, 14]
b_list = ['Python', 'Swift', 'Ruby', 'Go', 'Kotlin', 'Erlang']
b_list.sort()  # 对列表元素排序：默认按字符串包含的字符的编码大小比较
print(b_list)  # ['Erlang', 'Go', 'Kotlin', 'Python', 'Ruby', 'Swift']

# 指定key为len，指定使用len函数对集合元素生成比较的键，
# 也就是按字符串的长度比较大小
b_list.sort(key=len)
print(b_list)                       # ['Go', 'Ruby', 'Swift', 'Erlang', 'Kotlin', 'Python']
b_list.sort(key=len, reverse=True)  # 指定反向排序
print(b_list)                       # ['Erlang', 'Kotlin', 'Python', 'Swift', 'Ruby', 'Go']

stack = []
stack.append("fkit")          # 向栈中“入栈”3个元素
stack.append("crazyit")
stack.append("Charlie")
print(stack)         # ['fkit', 'crazyit', 'Charlie']
print(stack.pop())   # 第一次出栈：最后入栈的元素被移出栈
print(stack)         # ['fkit', 'crazyit']
print(stack.pop())   # 再次出栈
print(stack)         # ['fkit']

a_list = ['crazyit', 20, -1.2]
a_tuple = tuple(a_list)   # 将列表转换成元组
print(a_tuple)

a_range = range(1, 5)       # 使用range()函数创建区间（range）对象
print(a_range)             # range(1, 5)
b_tuple = tuple(a_range)  # 将区间转换成元组
print(b_tuple)           #[1, 2, 3, 4]

c_tuple = tuple(range(4, 20, 3))   # 创建区间时还指定步长
print(c_tuple)                    # [4, 7, 10, 13, 16, 19]

a_list = [2, 4, -3.4, 'crazyit', 23]
a_list[2] = 'fkit'        # 对第3个元素赋值
print(a_list)             # [2, 4, 'fkit', 'crazyit', 23]
a_list[-2] = 9527         # 对倒数第2个元素赋值
print(a_list)             # [2, 4, 'fkit', 9527, 23]

b_list = list(range(1, 5))
print(b_list)
b_list[1: 3] = ['a', 'b']   # 将第2个到第4个（不包含）元素赋值为新列表的元素
print(b_list)               # [1, 'a', 'b', 4]
b_list[1: 3] = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
print(b_list) # [1, 'a', 'b', 4]


# 将第3个到第3个（不包含）元素赋值为新列表的元素，就是插入
b_list[2: 2] = ['x', 'y']
print(b_list) # [1, 'a', 'x', 'y', 'b', 4]

# 将第3个到第6个（不包含）元素赋值为空列表，就是删除
b_list[2: 5] = []
print(b_list) # [1, 'a', 4]

# Python会自动将str分解成序列
b_list[1: 3] = 'Charlie'
print(b_list) # [1, 'C', 'h', 'a', 'r', 'l', 'i', 'e']

c_list = list(range(1, 10))
# 指定step为2，被赋值的元素有4个，因此用于赋值的列表也必须有4个元素
c_list[2: 9: 2] = ['a', 'b', 'c', 'd']
print(c_list) # [1, 2, 'a', 4, 'b', 6, 'c', 8, 'd']

d_list = list(range(1, 11))
d_list[2: 5: 1] = [80, 81, 82, 83]
d_list = list(range(1, 11))
##d_list[2: 5: 2] = [80, 81, 82, 83]
d_list[2: 5: 2] = [80, 81]
print('##6.end####################################################################################\n\n\n\n\n')




print('##7.start####################################################################################\n')
cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
print(cars['BMW'])
print(cars)   # {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
cars.clear()  # 清空cars所有key-value对
print(cars)   # {}

scores = {'语文': 89, '数学': 92, '英语': 93}
print(scores)
empty_dict = {}
print(empty_dict)  # 空的花括号代表空的dict
dict2 = {(20, 30):'good', 30:'bad'}   # 使用元组作为dict的key
print(dict2)



vegetables = [('celery', 1.58), ('brocoli', 1.29), ('lettuce', 2.19)]
# 创建包含3组key-value对的字典
dict3 = dict(vegetables)
print(dict3) # {'celery': 1.58, 'brocoli': 1.29, 'lettuce': 2.19}
cars = [['BMW', 8.5], ['BENS', 8.3], ['AUDI', 7.9]]
# 创建包含3组key-value对的字典
dict4 = dict(cars)
print(dict4) # {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
# 创建空的字典
dict5 = dict()
print(dict5) # {}
# 使用关键字参数来创建字典
dict6 = dict(spinach = 1.39, cabbage = 2.59)
print(dict6) # {'spinach': 1.39, 'cabbage': 2.59}


scores = {'语文': 89}
# 通过key访问value
print(scores['语文'])
# 对不存在的key赋值，就是增加key-value对
scores['数学'] = 93
scores[92] = 5.7
print(scores) # {'语文': 89, '数学': 93, 92: 5.7}
# 使用del语句删除key-value对
del scores['语文']
del scores['数学']
print(scores) # {92: 5.7}

cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
# 对存在的key-value对赋值，改变key-value对
cars['BENS'] = 4.3
cars['AUDI'] = 3.8
print(cars) # {'BMW': 8.5, 'BENS': 4.3, 'AUDI': 3.8}

# 判断cars是否包含名为'AUDI'的key
print('AUDI' in cars) # True
# 判断cars是否包含名为'PORSCHE'的key
print('PORSCHE' in cars) # False
print('LAMBORGHINI' not in cars) # True

print( dir(dict) )

print('##7.end####################################################################################\n\n\n\n\n')




print('##8.start####################################################################################\n')
# 字符串模板中使用key of dict
template1 = '书名是:%(name)s, 价格是:%(price)010.2f, 出版社是:%(publish)s'
book = {'name':'疯狂Python讲义', 'price': 88.9, 'publish': '电子社'}

# 使用字典为字符串模板中的key传入值
print(template1 % book)

book = {'name':'疯狂Kotlin讲义', 'price': 78.9, 'publish': '电子社'}
# 使用字典为字符串模板中的key传入值
print(template1 % book)

# 使用列表创建包含2个key的字典
a_dict = dict.fromkeys(['a', 'b'])
print(a_dict) # {'a': None, 'b': None}
# 使用元组创建包含2个key的字典
b_dict = dict.fromkeys((13, 17))
print(b_dict) # {13: None, 17: None}
# 使用元组创建包含2个key的字典，指定默认的value
c_dict = dict.fromkeys((13, 17), 'good')
print(c_dict) # {13: 'good', 17: 'good'}

cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
# 获取'BMW'对应的value
print(cars.get('BMW')) # 8.5
print(cars.get('PORSCHE')) # None
##print(cars['PORSCHE']) # KeyError

cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
# 获取字典所有的key-value对，返回一个dict_items对象
ims = cars.items()
print(type(ims)) # <class 'dict_items'>
# 将dict_items转换成列表
print(list(ims)) # [('BMW', 8.5), ('BENS', 8.3), ('AUDI', 7.9)]
# 访问第2个key-value对
print(list(ims)[1]) # ('BENS', 8.3)

# 获取字典所有的key，返回一个dict_keys对象
kys = cars.keys()
print(type(kys)) # <class 'dict_keys'>
# 将dict_keys转换成列表
print(list(kys)) # ['BMW', 'BENS', 'AUDI']
# 访问第2个key
print(list(kys)[1]) # 'BENS'

# 获取字典所有的value，返回一个dict_values对象
vals = cars.values()
# 将dict_values转换成列表
print(type(vals)) # [8.5, 8.3, 7.9]
# 访问第2个value
print(list(vals)[1]) # 8.3

print('##8.end####################################################################################\n\n\n\n\n')




print('##9.start####################################################################################\n')
cars = {'AUDI': 7.9, 'BENS': 8.3, 'BMW': 8.5}
print(cars)
# 弹出字典底层存储的最后一个key-value对
print(cars.popitem()) # ('AUDI', 7.9)
print(cars) # {'BMW': 8.5, 'BENS': 8.3}
# 将弹出项的key赋值给k、value赋值给v
k, v = cars.popitem()
print(k, v) # BENS 8.3


a = list(range(1,11))
print(a)
print(a[:])

cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
print(cars)
print(cars.pop('AUDI')) # 7.9
print(cars) # {'BMW': 8.5, 'BENS': 8.3}

cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
print(cars)
# 设置默认值，该key在dict中不存在，新增key-value对
print(cars.setdefault('PORSCHE', 9.2)) # 9.2
print(cars)
# 设置默认值，该key在dict中存在，不会修改dict内容
print(cars.setdefault('BMW', 3.4)) # 8.5
print(cars)

cars = {'BMW': 8.5, 'BENS': 8.3, 'AUDI': 7.9}
print(cars)
cars.update({'BMW':4.5, 'PORSCHE': 9.3})
print(cars)
print('##9.end####################################################################################\n\n\n\n\n')











