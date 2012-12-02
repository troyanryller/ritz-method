#! /usr/bin/env ruby
#encoding: UTF-8
require "../lib/main.rb"

describe Cmmffirst do
	it "h_i" do
		test = Cmmffirst.new(-1.0, 1.0, 5)
		test.h_i(3).should >= ( test.h - 0.00000001 )
	end
end