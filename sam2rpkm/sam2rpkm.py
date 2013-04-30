#!/usr/bin/python

################################################################################
#                                                                              #
# This script makes RPKM calculation in a mapping file (format sam) using a    #
# reference genome (genkbank format)                                           #
#                                                                              #
################################################################################

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
##any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = 'Jose Ricardo Bustos Molina'
__copyright__ = 'Copyright 2012, Jose Ricardo Bustos Molina, Corpogen'
__license__ = 'GPLv3 or later'
__version__ = '0.1.3'
__email__ = 'jrbustosm@gmail.com'
__status__ = 'beta'

import struct
import sys
import os

#Esta funcion se usa para comparar la posicion de un read respecto a un gen
#-1 si esta antes, 0 si esta adentro y 1 si esta despues
#para que un read se considere adentro de un gen al menos el 50% del read
#debe estar adentro del gen
#Para este algoritmo asumi que es imposible que el tamanyo de un read sea superior
#al tamnyo de cualquier gen
def comparar(a,b):
  if a[2] < b[1]:
    #Si el extremo derecho del read es inferior al extremo izquierdo del gen esta antes
    return -1
  if a[1] > b[2]:
    #Si el extremo izquierdo del read es superior al extremo derecho del gen esta despues
    return 1
  if a[1] >= b[1] and a[2] <= b[2]:
    #Si los dos extremos del read estan adentro del gen
    return 0
  if a[1] < b[1]:
    if a[2]-b[1]>=(a[2]-a[1])*50/100:
      #si el extremo izquierdo del read esta por fuera pero el extremo derecho esta por dentro
      #al menos el 50% del read debe estar dentro del gen
      return 0
    else:
      return -1
  else:
    if b[2]-a[1]>=(a[2]-a[1])*50/100:
      #si el extremo derecho del read esta por fuera pero el extremo izquierdo esta por dentro
      #al menos el 50% del read debe estar dentro del gen
      return 0
    else:
      return 1

def comparar2(a,b):
  if a[2] < b[1]:
    #Si el extremo derecho del read es inferior al extremo izquierdo del gen esta antes
    return -1
  if a[1] > b[2]:
    #Si el extremo izquierdo del read es superior al extremo derecho del gen esta despues
    return 1
  if a[1] >= b[1] and a[2] <= b[2]:
    #Si los dos extremos del read estan adentro del gen
    return 0
  if b[1] >= a[1] and b[2] <= a[2]:
    #el gen esta adentro del read, considero que hace parte del gen
    return 0
  if a[1] < b[1]:
    if a[2]-b[1]>=(a[2]-a[1])*50/100:
      #si el extremo izquierdo del read esta por fuera pero el extremo derecho esta por dentro
      #al menos el 50% del read debe estar dentro del gen
      return 0
    else:
      return -1
  else:
    if b[2]-a[1]>=(a[2]-a[1])*50/100:
      #si el extremo derecho del read esta por fuera pero el extremo izquierdo esta por dentro
      #al menos el 50% del read debe estar dentro del gen
      return 0
    else:
      return 1

def read_config():
  '''It reads the configuration options from the command line arguments and it returns a dict with them.'''
  from optparse import OptionParser
  usage = "usage: %prog [options] <sam> <gb>\nusage: %prog -c -l <INT> -u <int> <sam>"
  parser = OptionParser(usage=usage)
  parser.add_option("-o", "--out", dest="out", default="", help="Archivo de Salida, por defecto salida estandar", metavar="FILE")
  parser.add_option("-r", "--exclude", dest="exc", default="", help="Lista de genes a excluir en el analisis", metavar="FILE")
  parser.add_option("-f", "--filter", dest="filter", default="CDS", help="Filtro de los datos, por defecto CDS, ejemplo: 'CDS,tRNA'", metavar="STRING")
  parser.add_option("-l", "--lower", dest="lower", default="", help="Limite inferior a contar", metavar="INT")
  parser.add_option("-u", "--upper", dest="upper", default="", help="Limite superior a contar", metavar="INT")
  parser.add_option("-c", "--count", action="store_true", dest="count", default=False, help="solo contar, si tiene habilitados lower y upper, solo cuenta esa region")
  parser.add_option("-F", "--Flu", dest="flu", default="", help="Archivos con limites inferior y superior a contar", metavar="INT")
  parser.add_option("-i", "--intergenic", action="store_true", dest="intergenic", default=False, help="imprimir zonas intergenicas")
  parser.add_option("-I", "--intergenic_strict", action="store_true", dest="intergenic_strict", default=False, help="imprimir zonas intergenicas estrictas")
  parser.add_option("-L", "--log", action="store_true", dest="log", default=False, help="muestra una lista separada por comas de las secuencias encontradas en los calculos")
  parser.add_option("-a", "--antisense", action="store_true", dest="antisense", default=False, help="imprimir zonas antisentido")
  parser.add_option("-B", "--bruta", action="store_true", dest="bruta", default=False, help="hacer fuerza bruta cuando flu esta activado")

  #we parse the cmd line
  (options, args) = parser.parse_args()

  #we put the result in a dict
  global config
  config = {}
  for property in dir(options):
    if property[0] == '_' or property in ('ensure_value', 'read_file', 'read_module'): continue
    config[property] = getattr(options, property)

  if config['lower'] != '' or config['upper'] != '' or config['flu']:
    if not config['count'] or len(args)!=1:
      sys.argv.append('-h')
      read_config()
      sys.exit()
  else:   
    if len(args)==2:
      config['gb'] = args[1]
    else:
      sys.argv.append('-h')
      read_config()
      sys.exit()

  config['sam'] = args[0] 
  return config, args

def main():
  argv = sys.argv
  config, args = read_config()
  if len(args) == 0:
    sys.argv.append('-h')
    read_config()
    sys.exit()
  p = open(config['sam'])
  ls0 = []
  ls16 = []
  for l in p:
    if l[0] == '@':
      continue
    sp = l.split("\t")
    if(sp[1]=='0'):
      #Secuencias que mapearon por la hebra principal
      ls0.append((sp[0].split()[0], int(sp[3])-1, len(sp[9])+int(sp[3])-2))
    elif sp[1]=='16':
      #Secuencias que mapearon por la hebra complementaria
      ls16.append((sp[0].split()[0], int(sp[3])-1, len(sp[9])+int(sp[3])-2))

  p.close()
 
  #Si lo unico que se quiere es contar el numero de reads de una zona en particular usando un archivo de zonas
  if config['count'] and config['flu']!='':
    point = open(config['flu'])
    lines = point.readlines()
    zs = [('zone', int(l.split()[0]), int(l.split()[1])) for l in lines]

    if config['bruta']:
      for z in zs:
        s1 = sum([1 for v in ls0 if comparar2(v, z)==0])
        s2 = sum([1 for v in ls16 if comparar2(v, z)==0])
        print str(z[1]) + "\t" + str(z[2]) + "\t" + str(s1) + "\t" + str(s2)

      sys.exit()

    from operator import itemgetter
    zs.sort(key=itemgetter(1))

    def ordenar(a,b):
      c = (float(a[1]+a[2])/2)-(float(b[1]+b[2])/2)
      if c < 0:
        return -1
      if c > 0:
        return 1
      return 0

    ls0.sort(cmp=ordenar)
    ls16.sort(cmp=ordenar)
    i=0
    j=0
    m0 = dict([(z,0) for z in zs])
    while i<len(zs):
      c=comparar2(ls0[j],zs[i])
      if c==-1:
        j+=1
      elif c==0:
        j+=1
        m0[zs[i]]+=1
      else:
        i+=1

      if j>=len(ls0):
        break;

    i=0
    j=0
    m16 = dict([(z,0) for z in zs])
    while i<len(zs):
      c=comparar2(ls16[j],zs[i])
      if c==-1:
        j+=1
      elif c==0:
        j+=1
        m16[zs[i]]+=1
      else:
        i+=1

      if j>=len(ls16):
        break;

    for z in zs:
      print str(z[1]) + "\t" + str(z[2]) + "\t" + str(m0[z]) + "\t" + str(m16[z])

    sys.exit()

  #Si lo unico que se quiere es contar el numero de reads de una zona en particular
  if config['count'] and config['lower']!='' and config['upper']!='':
    z = ('zone', int(config['lower']), int(config['upper']))
    if not config['log']:
      s1 = sum([1 for v in ls0 if comparar2(v, z)==0])
      s2 = sum([1 for v in ls16 if comparar2(v, z)==0])
    else:
      s1=0
      s2=0
      log = []
      log2 = []
      for v in ls0:
        if comparar2(v, z)==0 :
          s1+=1
          log.append(v[0])
      for v in ls16:
        if comparar2(v, z)==0 :
          s2+=1
          log2.append(v[0])
    total = len(ls0) + len(ls16)
    print "Hebra principal: " + str(s1)
    print "Hebra complementaria: " + str(s2)
    print "total de reads: " + str(total )
    print "porcentaje: " + str((s1+s2) * 100.0 / total) + "%"
    if config['log']:
      print "Reads hebra principal:"
      print ','.join(log)
      print "Reads hebra complementaria:"
      print ','.join(log2)
    sys.exit()

  from operator import itemgetter
  #Ordenamos con base a la posicion de alineamiento
  ls0=sorted(ls0, key=itemgetter(1))
  ls16=sorted(ls16, key=itemgetter(1))

  from Bio import SeqIO
  info = SeqIO.parse(open(config['gb']), 'gb').next()
  fs = [f for f in info.features if f.type in config['filter']]
  #Filtramos genes ribosomales
  if config['exc']!='':
    p = open(config['exc'])
    rbs = [l.strip() for l in p.readlines()]
    fs = [f for f in fs if not f.qualifiers['locus_tag'][0] in rbs]

  if not config['antisense']:
    #genes en la hebra principal
    fs1 = [f for f in fs if f.strand==1]
    #genes en la hera complementaria
    fs_1 = [f for f in fs if f.strand==-1]
  else:
    #cambiamos las dos hebras de manera predeterminada
    fs1 = [f for f in fs if f.strand==-1]
    fs_1 = [f for f in fs if f.strand==1]
  del fs

  #Buscamos las posiciones de cada gen
  fs1 = [(f.qualifiers['locus_tag'][0], f.location.start.position, f.location.end.position-1) for f in fs1]
  fs_1 = [(f.qualifiers['locus_tag'][0], f.location.start.position, f.location.end.position-1) for f in fs_1]

  dfs1 = dict([(d[0], d) for d in fs1])
  dfs_1 = dict([(d[0], d) for d in fs_1])

  if config['intergenic_strict']:
    mapp = [True for i in range(len(info))]
    for v in fs1+fs_1:
      for i in range(v[1],v[2]+1):
        mapp[i]=False

    flag = True
    lst = []
    l=0
    u=0
    for i in range(len(mapp)):
      if mapp[i] == flag:
        if flag:
          l = i
        else:
          u = i-1
          lst.append( ('int_strict' + str(len(lst)), l+1, u+1) )
        flag = not flag
    for z in lst:
      if not config['log']:
        s1 = sum([1 for v in ls0 if comparar(v, z)==0])
        s2 = sum([1 for v in ls16 if comparar(v, z)==0])
        print z[0] + "\t" + str(z[1]) + "-" + str(z[2]) + "\t" + str(s1) + "\t" +  str(s2) + "\t" + str(s1+s2)
      else:
        s1=0
        s2=0
        log = []
        log2 = []
        for v in ls0:
          if comparar(v, z)==0 :
            s1+=1
            log.append(v[0])
        for v in ls16:
          if comparar(v, z)==0 :
            s2+=1
            log2.append(v[0])
        print z[0] + "\t" + str(z[1]) + "-" + str(z[2]) + "\t" + str(s1) + "\t" +  str(s2) + "\t" + str(s1+s2) + "\t" + ",".join(log) + "\t" + ",".join(log2)
    sys.exit()

  j=0 #Esta variable lleva la cuenta en que gen vamos
  mapGen = {} #Almacena los read encontrados para cada gen
  mapInt = {} #Almacena los read encontrados para cada gen en una zona intergenica
  for read in ls0:  
    if(len(fs1)==0):
      break

    while j<len(fs1)-1 and comparar(read, fs1[j]) > 0:
      j+=1

    if(comparar(read, fs1[j])==0):
      if not (fs1[j][0] in mapGen):
        mapGen[fs1[j][0]] = []
      mapGen[fs1[j][0]].append(read)
    elif(comparar(read, fs1[j-1])==0):
      #caso especial: a veces hay un problema con el ordenamiento y queda una secuencia corta
      #que pertenece al gen anterior marcada para el siguiente gen
      if not (fs1[j-1][0] in mapGen):
        mapGen[fs1[j-1][0]] = []
      mapGen[fs1[j-1][0]].append(read)
    else:
      #encontramos un nombre apropiado para la zona intergenica
      if j==0:
        intname = "1-" + str(fs1[0][1])
      elif j==len(fs1)-1:
        intname = str(fs1[j][2]+2) + "-" + str(len(info))
      else:
        intname = str(fs1[j-1][2]+2) + "-" + str(fs1[j][1])
      #Agregamos a mapa de Zonas intergenicas
      if not (intname in mapInt):
        mapInt[intname] = []
      mapInt[intname].append(read)

  j=0 #Esta variable lleva la cuenta en que gen vamos
  mapGen2 = {} #Almacena los read encontrados para cada gen de la hebra complemntaria
  mapInt2 = {} #Almacena los read encontrados para cada gen en una zona intergenica
  for read in ls16:  
    if(len(fs_1)==0):
      break

    while j<len(fs_1)-1 and comparar(read, fs_1[j]) > 0:
      j+=1

    if(comparar(read, fs_1[j])==0):
      #Si el read esta adentro de un gen lo anyadimos
      if not (fs_1[j][0] in mapGen2):
        mapGen2[fs_1[j][0]] = []
      mapGen2[fs_1[j][0]].append(read)
    elif(comparar(read, fs_1[j-1])==0):
      #caso especial: a veces hay un problema con el ordenamiento y queda una secuencia corta
      #que pertenece al gen anterior marcada para el siguiente gen
      if not (fs_1[j-1][0] in mapGen2):
        mapGen2[fs_1[j-1][0]] = []
      mapGen2[fs_1[j-1][0]].append(read)
    else:
      #encontramos un nombre apropiado para la zona intergenica
      if j==0:
        intname = "1-" + str(fs_1[0][1])
      elif j==len(fs_1)-1:
        intname = str(fs_1[j][2]+2) + "-" + str(len(info))
      else:
        intname = str(fs_1[j-1][2]+2) + "-" + str(fs_1[j][1])
      #Agregamos a mapa de Zonas intergenicas
      if not (intname in mapInt2):
        mapInt2[intname] = []
      mapInt2[intname].append(read)

  if not config['intergenic']:
    rpkm = {}
    for k,v in mapGen.iteritems():
      if not config['count']:
        rpkm[k] = len(v)/(((len(ls16)+len(ls0)) / 1000000.0)*((dfs1[k][2]-dfs1[k][1])/1000.0))
      else:
        rpkm[k] = len(v)

    rpkm2 = {}
    for k,v in mapGen2.iteritems():
      if not config['count']:
        rpkm2[k] = len(v)/(((len(ls16)+len(ls0)) / 1000000.0)*((dfs_1[k][2]-dfs_1[k][1])/1000.0))
      else:
        rpkm2[k] = len(v)

    for k in dfs1:
      if k in rpkm:
        aux = ""
        if config['log']:
          aux += "\t" + ",".join([v[0] for v in mapGen[k]])
        print k + "\t" + str(rpkm[k]) + aux
      else:
        print k + "\t0"

    for k in dfs_1:
      if k in rpkm2:
        aux = ""
        if config['log']:
          aux += "\t" + ",".join([v[0] for v in mapGen2[k]])
        print k + "\t" + str(rpkm2[k]) + aux
      else:
        print k + "\t0"
  else:
    rpkm = {}
    for k,v in mapInt.iteritems():
      l = k.split("-")
      rpkm[k] = len(v)/(((len(ls16)+len(ls0)) / 1000000.0)*((float(l[1])-float(l[0])+2)/1000.0))
      aux = ""
      if config['log']:
        aux += "\t" + ",".join([d[0] for d in v])
      print k + "\t" + str(rpkm[k]) + "\t" + str(len(v)) + "\t" + "1" + aux

    rpkm2 = {}
    for k,v in mapInt2.iteritems():
      l = k.split("-")
      rpkm2[k] = len(v)/(((len(ls16)+len(ls0)) / 1000000.0)*((float(l[1])-float(l[0])+2)/1000.0))
      aux = ""
      if config['log']:
        aux += "\t" + ",".join([d[0] for d in v])
      print k + "\t" + str(rpkm2[k]) + "\t" + str(len(v)) + "\t" + "-1" + aux

if __name__ == "__main__":
  sys.exit(main())

