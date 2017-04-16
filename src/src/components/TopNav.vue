<template>
  <el-row class="top-nav">
    <el-menu default-active="1" class="el-menu-demo" mode="horizontal" @select="handleSelect">
      <el-col :span="22" :offset="1">
        <el-menu-item index="1" v-if="'logo' in top" >
          <img class="nsp-logo" :src="'static/'+top.logo" />
        </el-menu-item>
        <el-menu-item index="1">
          <router-link to="/home">Home</router-link>
        </el-menu-item>
        <el-menu-item index="1">
          <router-link to="/doc">Documentation</router-link>
        </el-menu-item>
        <!--<el-menu-item index="1">
          <router-link to="/api">API</router-link>
        </el-menu-item>-->
        <el-menu-item index="1">
          <router-link to="/download">Download</router-link>
        </el-menu-item>
        <el-menu-item index="1">
          <el-menu-item index="3"><a href="mailto:wj_hust08@hust.edu.cn?cc=yxiao@hust.edu.cn">Contact</a></el-menu-item>
      </el-col>
    </el-menu>
    <div class="line"></div>
  </el-row>
</template>

<script>
  export default {
    name: 'nav',

    data() {
      return {
        top: []
      }
    },

    created() {
      this.fetch_data()
    },

    methods: {
      fetch_data() {
        var v = this
        this.$http.get('static/top.json').then((response) => {
          v.top = response.body
        }, (response) => {
        })
      },

      handleSelect(key, keyPath) {
        // console.log(key, keyPath)
      }
    }
  }
</script>

<style scoped>
  a {
    color: rgb(71, 86, 105);
    text-decoration: none;
  }

  .top-nav {
    border-bottom: 1px solid #f0f0f0;
  }

  .el-menu {
    background-color: #fff;
  }

  .nsp-logo {
    vertical-align: middle;
    height: 50px;
  }

</style>
