<template>
  <el-row>
    <el-col :span="16" :offset="4">
      <el-card class="download">
        <ul>
          <li><a href="http://biophy.hust.edu.cn/download/nsp.tar.gz">Source code</a></li>
          <li><a href="http://biophy.hust.edu.cn/download/nsp-lib.tar.gz">Templates library</a></li>
        </ul>
      </el-card>
    </el-col>
  </el-row>
</template>

<script>
  import { bus } from '../bus.js'

  export default {
  }
</script>

<style>
  .download {
    margin: 10px;
    padding: 20px;
  }
</style>