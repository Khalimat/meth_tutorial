import numpy as np
from scipy.signal import argrelextrema

'''Создаем словарь, где в качестве key служат названия хромосом, а значений — трехмерный массив 
(val[0]- позиции, val[1] - значения градиента метилирования, val[2]- значения градиента стандартного отклонения, 
val[3] - промежуток)'''
results = dict((key, [[], [], [], []]) for key in ["chr" + str(i) for i in range(1, 23)])

'''Создаем словарь, где в качестве ключей служат названия хромосом, а значений — трехмерный массив. 
Мы будем записывать в значения трехмерный маcсив, где val[0]- градиент среднего метилирования, val[1]-градиент среднего 
стандатного отклонения метилирования, val[2]- позицию на хромосоме
'''
gradient = dict((key, [[], [], []]) for key in ["chr" + str(i) for i in range(1, 23)])


'''Создаем словарь, где в качестве ключей служат названия хромосом, а значений — трехмерный массив, 
где val[0]- average_meth, val[1]-variance, val[2]-start_position
'''
genome_coordinates = dict((key, [[], [], []]) for key in ["chr" + str(i) for i in range(1, 23)])


'''Функция для пересечения двух списков, выдает список общих элементов'''
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

'''zip позволяет пройтись одновременно по нескольким итерируемым объектам. Ходим сразу по элементам двух файлов.'''
with open("/home/khali/IOGen/MetAverage200.txt", "r") as foo, open("/home/khali/IOGen/MetStdev200.txt", "r") as bar:
    for f, b in zip(foo, bar):
        '''Читаем построчно файлы и сохраняем номер хромосомы, позицию старта и значение среднего метилирования
               из файла MetAverage200'''
        line_av_chr, line_av_start, line_av_value = f.strip().split()[0:2] + f.strip().split()[3:]

        '''Cохраняем номер хромосомы, позицию старта и значение var метилирования 
        из файла MetStdev200.txt'''
        line_var_chr, line_var_start, line_var_value = b.strip().split()[0:2] + b.strip().split()[3:]


        '''Проверяем совпадают ли названия хромосомы и позиции старта в двух файлах, они должны совпадать, но 
        всегда лучше все перепроверить. Если да, то добавляем в genome_coordinates[chrN] последовательно average, 
        variance и позицию старта'''
        if line_av_chr == line_var_chr and line_av_start == line_var_start:
            genome_coordinates[line_av_chr][0].append(line_av_value)
            genome_coordinates[line_var_chr][1].append(line_var_value)
            genome_coordinates[line_var_chr][2].append(line_var_start)

    '''Последовательно перебираем значения в словаре genome_coordinates, создаем np-массивы из массивов со 
        значениями average и var, вычисляем для каждого элемента этих массивов градиент. Записываем последовательно
        массивы c grad_average, grad_var и позициями в словарь gradient, ключами служат названия хромосом'''
    for key, val in genome_coordinates.items():
        np_average = np.array(val[0], dtype=float)
        np_gradient_average = np.gradient(np_average)
        np_var = np.array(val[1], dtype=float)
        np_gradient_var = np.gradient(np_var)
        np_pos = np.array(val[2], dtype=int)
        gradient[key][0].append(np_gradient_average)
        gradient[key][1].append(np_gradient_var)
        gradient[key][2].append(np_pos)

    '''Находим минимумы grad в окрестности +/- 20 b.p. в variance и average c помощью функции argrelextrema (возвращает индексы), 
    пересекаем массивы с полученными индексами минимумов average и variance. Добавляем в results позиции полученные в результате пересечения,
     и значения градиента average и var в этой позиции. Строчки 74-75 добавляет промежуток +-20 СpG при условии, что это не граничный участок'''
for key, val in gradient.items():
    average = val[0][0]
    variance = val[1][0]
    min_average = argrelextrema(average, np.less, order=20) 
    min_variance = argrelextrema(variance, np.less, order=20)
    indexes = intersection(min_average[0], min_variance[0])
    for i in indexes:
        results[key][0].append(gradient[key][2][0][i])
        results[key][1].append(gradient[key][0][0][i])
        results[key][2].append(gradient[key][1][0][i])
        if i > 20 and len(gradient[key][2][0]) - i > 20:
            z = str(key) + ":" + str(gradient[key][2][0][i - 20]) + "-" + str(gradient[key][2][0][i + 20])
            results[key][3].append(z)

'''Записываем в obtained_sites_2.txt номера хромосом, позиции, min(grad_var), min(grad_average) и целевой промежуток'''
with open("/home/khali/IOGen/target_sites_grad.txt", "w") as target:
    print("chr_ID", "pos", "grad_meth_av", "grad_methy_var",
          sep='\t', file=target)

    for key, val in results.items():
        for j in range(len(val[0])):
            print(key, val[0][j], np.around(val[1][j], decimals=9), np.around(val[2][j], decimals=9), val[3][j],
                  sep='\t', file=target)
