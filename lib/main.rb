#! /usr/bin/env ruby
#encoding: UTF-8

module Math
	#Getting integral by the method of rectangles
	 #a - lower limit of integration
	 #b - upper limit of integration
	 #n - number of iterations
	def integral a, b, n = 10000	
		h = (b - a).to_f / n.to_f
		x_i = []
		(n+1).times { |i| x_i << (h*i + a).to_f }
		sum = 0.0
		(n+1).times { |i| sum += yield(x_i[i] - h/2) }
		h * sum
	end
end

class Cmmffirst
	include Math
	attr_reader :h, :x_i

  def initialize(a,b,n)
    @a, @b, @n = a, b, n
    @h = ( @b.to_f - @a.to_f ) / @n.to_f
    @x_i = []
    (@a..@b).step(@h) { |x| @x_i << x }
    #puts @x_i[0]
    #(@n+1).times { |i| @point_x_i << @h*i }
  end

  def h_i i
    i+=1
    @x_i[i] - @x_i[i-1]
  end

  def p x
    #print "Enter function p(x): "
    #str = gets
    #eval str
    (x + 2.0).to_f / (3.0 * x + 4.0).to_f
  end

  def q x
    #print "Enter function q(x): "
    #str = gets
    #eval str
    1 + sin(x)
  end

  def p_derivative x
    #print "Enter derivative of function p(x): "
    #str = gets
    #eval str
    10 / (9 * x * x + 24 * x + 16)
  end

  def q_derivative x
    #print "Enter derivative of function p(x): "
    #str = gets
    #eval str
    cos(x)
  end

    def make_pfi(x, i)
    case x
      #when i == @n - 1
      #  0.0
    when @x_i[i-1] ... @x_i[i] 
          ( x - @x_i[i-1] ) / h_i(i).to_f
    when @x_i[i] ... @x_i[i+1] 
        - ( x - @x_i[i+1] ) / h_i(i+1).to_f
    
    else
      0.0
    end
  end

  def make_derivative_pfi(x, i)
    case x
    when @x_i[i-1] ... @x_i[i] 
          1.0 / h_i(i).to_f 
    when @x_i[i] ... @x_i[i+1] 
        - 1.0 / h_i(i+1).to_f
    else
      0.0
    end
  end

  def scalar_p_pfipfi i, j
    i+=1
    j+=1
    integral(@x_i[i-1], @x_i[i], 10000) { |x| ( p(x) * make_derivative_pfi(x, i) * make_derivative_pfi(x, j) )  }
  end

  def scalar_q_pfipfi i, j
    i+=1
    j+=1
    integral(@x_i[i-1], @x_i[i], 10000) { |x| ( q(x) * make_pfi(x, i) * make_pfi(x, j) )  }
    #@x_i[i-1]
  end

  
  def matrix_m_i i
    Matrix.columns([ [scalar_q_pfipfi(i-1, i-1), scalar_q_pfipfi(i-1, i)], 
                     [scalar_q_pfipfi(i, i-1), scalar_q_pfipfi(i, i)] ])
  end

  def matrix_k_i i
    Matrix.columns([ [scalar_p_pfipfi(i-1, i-1), scalar_p_pfipfi(i-1, i)], 
                     [scalar_p_pfipfi(i, i-1), scalar_p_pfipfi(i, i)] ])
    #scalar_p_pfipfi(i-1, i)
  end

  def matrix_a i
    #rescue "err" if ((i == 0) || (i == @n-1)) 
    rez = []
    tmp = []
    tmp << ( scalar_q_pfipfi(i-1, i-1) + scalar_p_pfipfi(i-1, i-1) )
    tmp << ( scalar_q_pfipfi(i-1, i) + scalar_p_pfipfi(i-1, i) )
    #p tmp
    rez << tmp
    tmp = []
    tmp << ( scalar_q_pfipfi(i, i-1) + scalar_p_pfipfi(i, i-1) )
    tmp << ( scalar_q_pfipfi(i, i) + scalar_p_pfipfi(i, i) )
    #p tmp
    rez << tmp
    rez
  end

  def matrix_m_i_analit i
    raise Exception if i < 1
    element1 = ( 1.0 / ( h_i( i ) * h_i( i ) ) ) *
           ( 0.5 * ( @x_i[i] * @x_i[i] - @x_i[i-1] * @x_i[i-1] ) - 
             @x_i[i] * cos(@x_i[i]) + @x_i[i-1] * cos(@x_i[i-1]) - 
             sin(@x_i[i]) + sin(@x_i[i-1]) - @x_i[i]*@x_i[i] * ( @x_i[i] - @x_i[i-1] ) + 
             @x_i[i] * @x_i[i] * cos(@x_i[i]) - @x_i[i] * @x_i[i] * cos(@x_i[i-1]) )
    element2 = ( ( @x_i[i] * @x_i[i] * @x_i[i] ) / 3.0 - ( @x_i[i-1] * @x_i[i-1] * @x_i[i-1] ) / 3.0 -
           @x_i[i] * @x_i[i] * cos(@x_i[i]) + @x_i[i-1] * @x_i[i-1] * cos(@x_i[i-1]) + 
           2.0 * @x_i[i] * sin(@x_i[i]) - 2.0 * @x_i[i-1] * sin(@x_i[i-1]) + 
           2.0 * @x_i[i] * cos(@x_i[i]) - 2.0 * @x_i[i-1] * cos(@x_i[i-1]) - 
           ( ( @x_i[i-1] + @x_i[i] ) / 2.0 ) * ( @x_i[i]*@x_i[i] - @x_i[i-1]*@x_i[i-1] ) - 
           ( @x_i[i-1] + @x_i[i] ) * ( -@x_i[i] * cos(@x_i[i]) + @x_i[i-1] * cos(@x_i[i-1]) -
           sin(@x_i[i]) + sin(@x_i[i-1]) ) + @x_i[i] * @x_i[i-1] * ( @x_i[i] - @x_i[i-1] ) + 
           @x_i[i] * @x_i[i-1] * ( -cos(@x_i[i]) + cos(@x_i[i-1]) ) ) *( -1.0 )
    [ [element1, element2], 
      [element2, element1] ]
  end

  def matrix_k_i_analit i
    raise Exception if i < 1
    element = ( 1.0 / ( h_i( i ) * h_i( i ) ) ) *
           ( ( 1.0 / 3.0 ) * ( @x_i[i] - @x_i[i] ) + ( 2.0 / 9.0 ) * 
            ( log(3.0 * @x_i[i] + 4.0) - log(3.0 * @x_i[i-1] + 4.0) ) ) 
    [ [element, -element], 
      [-element, element] ]
  end

 def matrix_a_i_analit i
   #@n.times{ |i| ( matrix_k_i_analit(i+1) + matrix_m_i_analit(i+1) ) }
   Matrix.columns( matrix_k_i_analit(i) ) + Matrix.columns( matrix_m_i_analit(i) )
 end

 def matrix_a_analit 
  m = Matrix.zero(@n).to_a
  #m1 = matrix_a_i_analit(1).to_a
  @n.times do |i|  
    m2 = matrix_a_i_analit(i+1).to_a
    if i==0
      m[i][i] += m2[1][1]
    else
      m[i-1][i-1] += m2[0][0]
      m[i-1][i] += m2[0][1]
      m[i][i-1] += m2[1][0]
      m[i][i] += m2[1][1]
    end
  end
  me = Matrix[*m] 
 end
 
end

test = Cmmffirst.new(-1.0, 1.0, 5)
tmp1 = test.matrix_a_analit 
puts tmp1.row_vectors()
#p test.h
#p test.x_i[3] - test.x_i[2]
#p test.matrix_k_i_analit 1
#p test.matrix_k_i 1
#p test.matrix_m_i_analit 1
#p test.matrix_m_i 1
#p test.matrix_k_i 1
#p test.matrix_k_i 19
#5.times{ |i| p i }
#.times{ |i| p test.matrix_a(i) }
#p test.scalar_q_pfipfi(1, 1)
