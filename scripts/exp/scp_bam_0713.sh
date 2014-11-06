#!/bin/bash

while read s
do
	scp $URANUS:/rd1/user/liufl/software/apache-tomcat-6.0.37/ExpCall/bams/${s}.bam ./
	scp $URANUS:/rd1/user/liufl/software/apache-tomcat-6.0.37/ExpCall/bams/${s}.bam.bai ./
done
