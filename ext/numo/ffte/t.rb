require "numo/narray"
require "./ffte"

if true
#p n = 2**24
#p n = 12345
p n = 2**20
p a = Numo::DComplex.new(3,n).seq
p x = Numo::FFTE.zfft1d(a,1)
p y = Numo::FFTE.zfft1d(x,-1)
p (a-y).abs.max

p nx = 2**10
p ny = 2**10
p a = Numo::DComplex.new(3,ny,nx).seq
p x = Numo::FFTE.zfft2d(a,1)
p y = Numo::FFTE.zfft2d(x,-1)
p (a-y).abs.max

p nx = 2**7
p ny = 2**7
p nz = 2**7
p a = Numo::DComplex.new(3,ny,nx,nz).seq
p x = Numo::FFTE.zfft3d(a,1)
p y = Numo::FFTE.zfft3d(x,-1)
p (a-y).abs.max
end

p nx = 2**10/2+1
p ny = 2**10
p a = Numo::DComplex.new(3,ny,nx).seq
p x = Numo::FFTE.zdfft2d(a)
p y = Numo::FFTE.dzfft2d(x)
p (a-y).abs.max

if true
p nx = 2**7/2+1
p ny = 2**7
p nz = 2**7
p a = Numo::DComplex.new(3,nz,ny,nx).seq
p x = Numo::FFTE.zdfft3d(a)
p y = Numo::FFTE.dzfft3d(x)
p (a-y).abs.max
end
