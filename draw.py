# 随机抽签的程序，从一个字典中随机抽取一个名字，然后打印出来。
# 字典的内容来自students.txt文件，文件中每行第一列是名字，第二列是学号。
# 例如：
# 张三 20130001
# 李四 20130002
# 王老五 20130003

import random
import chardet

# 输出欢迎信息, 居中对齐
print('\n')
print('欢迎使用能源与动力工程系随机抽签程序v1.0 by: 火兴辉'.center(80))
print('\n')

file_name = 'students.txt'
# 读取文件内容
with open(file_name, 'rb') as f:
    data = f.read()
    file_encoding = chardet.detect(data)['encoding']

try:
    with open(file_name, 'r', encoding=file_encoding) as read_file:
        lines = read_file.readlines()
except UnicodeDecodeError:
    encodings_to_try = ['utf-8', 'gbk', 'gb2312', 'gb18030', 'big5']
    for file_encoding in encodings_to_try:
        try:
            with open(file_name, 'r', encoding=file_encoding) as read_file:
                lines = read_file.readlines()
                break
        except UnicodeDecodeError:
            continue
    else:
        print('文件编码错误，无法解码！')
        exit()

print('students.txt文件编码为：', file_encoding)

students = {}
for line in lines:
    name, number = line.split()
    students[name] = number

# 显示目前学生数量
print('目前学生数量：', len(students))

# 输入抽签人的姓名
drawer_name = input('请输抽签人的名字：')

# 输入抽签数量
draw_number = int(input('请输入抽签数量：'))

# 随机抽取draw_number个名字，如果抽签数量大于文件中的名字数量，则输出错误信息
if draw_number > len(students):
    print('抽签数量大于文件中的名字数量！')
    exit()
draw_names = random.sample(list(students.keys()), draw_number)

# 输出抽签结果并询问是否确认
print('抽签结果：')
for name in draw_names:
    print(name, students[name])
confirm = input('是否确认？(Y/N)')
# 如果不确认则退出程序, 否则写入文件, 忽略大小写
if confirm.lower() != 'y':
    exit()

# 写入文件, 以追加的方式写入, 如果文件不存在则创建。首先写入抽签人的名字，换行，然后每行写入一个抽取的名字和学号。
with open('draw_result.txt', 'a', encoding=file_encoding) as write_file:
    write_file.write(drawer_name + '\n')
    for name in draw_names:
        write_file.write(name + ' ' + students[name] + '\n')
    write_file.write('\n')
print('抽签结果已保存到draw_result.txt文件中。')

# 从students.txt文件中删除已经抽取的名字
with open(file_name, 'w', encoding=file_encoding) as write_file:
    for name in students.keys():
        if name not in draw_names:
            write_file.write(name + ' ' + students[name] + '\n')

# 等待用户按回车键退出程序
input('按回车键退出程序...')